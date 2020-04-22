function [signal_mod,y,opt_val_sort] = decomposeGauss5Cou(ppg_signal,freq,varargin)
% Pulse Decomposition Analysis: Gauss5 algorithm (Couceiro restrictions)
%
% References:
% Ricardo Couceiro u.a. ÑAssessment of cardiovascular function from
% multi-Gaussian fitting of a finger photoplethysmogramì. In: Physiological
% Measurement 36.9 (2015), S. 1801ñ1825.
%
% This function decomposes a single beat of a ppg signal by using 5 gaussian
% functions as basis functions with the restrictions of Couceiro et al. The 
% optimization algorithm sorts the basis functions by their respective 
% maxima location. The outputs are the 5 optimized basis functions, the 
% optimized parameters of the kernels as well as the reconstructed signal 
% which is the sum of the 5 basis functions.
% 
% input: 
% 'ppg_signal':     detrended ppg signal (single beat from minimum to 
%                   minimum) that is to be decomposed
%
% 'freq':           sampling frequency
%
% 'varargin':       optional arguments:
%                       'NoOpt' - logical, that determines whether (false -
%                       default) optimization is executed or the
%                       optimization is skipped (true)
%
%                       'InitialValues' - vector, that contains initial
%                       values for the kernels; if no vector is passed, the
%                       default initial values are used
%
% output:
% 'signal_mod':     reconstructed signal (sum of basis functions)
%
% 'y':              cell array of size 5x1 which contains the 5 basis
%                   functions (y1, y2, y3, y4, y5) each with a size of 
%                   1xlength(ppg_signal) sorted by the occurence of their
%                   maxima
%
% 'opt_values':     optimized values of parameters that determine the basis
%                   functions

%% check input arguments
% both signal and sampling frequency are needed
if(nargin<2)
    errordlg('Too few arguments','Input Error','modal');
    return;
elseif(nargin>6)
    errordlg('Too many arguments','Input Error','modal');
    return;
end
% ppg_signal needs to be an array with more than one element
if(~(isvector(ppg_signal) && ~(isscalar(ppg_signal)) && ~(ischar(ppg_signal)) && ~isstring(ppg_signal)))
    errordlg('PPG signal needs to be an array with more than one element',...
        'Input Error','modal');
    return;
end
% freq needs to be a scalar
if(~(isscalar(freq) && ~(ischar(freq)) && ~isstring(freq)))
    errordlg('Sampling frequency needs to be a scalar value',...
        'Input Error','modal');
    return;
end
% turn ppg_signal into 1xn array
if(size(ppg_signal,1)>size(ppg_signal,2))
    ppg_signal = ppg_signal';
end
% check optional arguments
okargs = {'InitialValues','NoOpt'} ;% allowed argument specifiers
% 'InitialValues' + vector - overwrite inital values stated in this
% function
% 'NoOpt' + logical - if true, optimization is skipped and model is created
% with inital values
x0 = []; % initial values
noOpt = []; % optimization flag
i=1;
while(1)
    if(mod(size(varargin,2),2)~=0) % check if number of variable arguments is even
        errordlg('Uneven number of optional arguments not supported',...
            'Input Error','modal');
        return;
    end
    
    if(i>size(varargin,2)) % leave while loop if there are no more new variable arguments
        break;
    end
    
    if(ischar(varargin{i})) % evaluate key word of next argument
        if(ismember(varargin{i},okargs))
            actualArgument=varargin{i};
        else
            errordlg(['Specified Option ''' varargin{i} ''' is not a valid option'],...
                'Input Error','modal');
            return;
        end
        
    else
        errordlg(['Argument ' num2str(i+2) ' is not a valid option. A keyword needs to be inserted here.'],...
            'Input Error','modal');
        return;
    end
    
    switch actualArgument % set values for the actual arguments
        case 'InitialValues'
            % turn ppg_signal into 1xn array
            if(size(varargin{i+1},1)>size(varargin{i+1},2))
                varargin{i+1} = varargin{i+1}';
            end
            if(~(isvector(varargin{i+1}) && (size(varargin{i+1},2)==9)))
                errordlg('InitialValues needs to be an array with a length of 9',...
                    'Input Error','modal');
                return;
            end
            x0=varargin{i+1};
        case 'NoOpt'
            if(~islogical(varargin{i+1}))
                errordlg('Optimization flag needs to be logical',...
                    'Input Error','modal');
                return;
            end
            noOpt=varargin{i+1};
    end
    i=i+2;
end

%% produce time axis for input beat
t_ppg= 0:1/freq:(length(ppg_signal)-1)/freq;

%% seperation into sytole and diastole
% finding dicrotic notch
deriv_2 = deriv2Couceiro(ppg_signal);

%do moving average filter to smoothen second derivative
deriv_2 = movmean(deriv_2,50);
%specify min peak prominence
minProm = (max(deriv_2)-min(deriv_2))/10;

sample_vec = 1:length(deriv_2);
interval_beg = 0.2*freq;
interval_end = 0.4*freq;
% find reference point for dicrotic notch
try
    [max_amp,max_loc] = findpeaks(deriv_2(interval_beg:interval_end),'MinPeakProminence',minProm);
    max_loc = max_loc+interval_beg;
    if(length(max_amp)>1)
        max_loc = max_loc(max_amp==max(max_amp));
        if(length(max_loc)>1) % if there is more than one peak of same height
            max_loc = max_loc(1); % take first peak
        end
        max_amp = max(max_amp);
    end
    if(isempty(max_loc))
        % if there are no peaks in temporal window, take maximum of signal
        % in defined interval
        max_amp = max(deriv_2(interval_beg:interval_end));
        max_loc = sample_vec(deriv_2 == max_amp);
        max_loc(max_loc>interval_end) = [];
        max_loc(max_loc<interval_beg) = [];
    end
catch
    % if there are no peaks in temporal window, take maximum of signal
    % in defined interval
    max_amp = max(deriv_2(interval_beg:interval_end));
    max_loc = sample_vec(deriv_2 == max_amp);
    max_loc(max_loc>interval_end) = [];
    max_loc(max_loc<interval_beg) = [];
end
dic_ref = max_loc;
% find boundaries
deriv_2_sign = sign(deriv_2);
zero_cross = sample_vec([0 diff(deriv_2_sign)] ~= 0);
zero_cross_before = zero_cross(zero_cross<max_loc);
zero_cross_before = max(zero_cross_before);
if(zero_cross_before<interval_beg)
    zero_cross_before = [];
else
    bound_before = zero_cross_before;
end
zero_cross_after = zero_cross(zero_cross>max_loc);
zero_cross_after = max(zero_cross_after);
if(zero_cross_after>interval_end)
    zero_cross_after = [];
else
    bound_after = zero_cross_after;
end
if(isempty(zero_cross_before))% needs refinement
   deriv_4 = deriv4Couceiro(ppg_signal);
   deriv_4_sign = sign(deriv_4);
   zero_cross4 = sample_vec([0 diff(deriv_4_sign)] ~= 0);
   zero_cross4_before = zero_cross4(zero_cross4<max_loc);
   zero_cross4_before = max(zero_cross4_before); 
   bound_before = zero_cross4_before;
end
if(isempty(zero_cross_after))
   deriv_4 = deriv4Couceiro(ppg_signal);
   deriv_4_sign = sign(deriv_4);
   zero_cross4 = sample_vec([0 diff(deriv_4_sign)] ~= 0);
   zero_cross4_after = zero_cross4(zero_cross4>max_loc);
   zero_cross4_after = min(zero_cross4_after); 
   bound_after = zero_cross4_after;
end

%% search for important points
[peaks_pos_amp,peaks_pos_loc] = findpeaks(deriv_2,'MinPeakProminence',minProm);
[peaks_neg_amp,peaks_neg_loc] = findpeaks(-deriv_2,'MinPeakProminence',minProm);
peaks_pos_loc_sys = peaks_pos_loc(peaks_pos_loc<=bound_before);
peaks_neg_loc_sys = peaks_neg_loc(peaks_neg_loc<=bound_before);
peaks_pos_loc_dias = peaks_pos_loc(peaks_pos_loc>=bound_after);
peaks_neg_loc_dias = peaks_neg_loc(peaks_neg_loc>=bound_after);
pos_a = peaks_pos_loc_sys(1);
pos_b = peaks_neg_loc_sys(1);
if(length(peaks_pos_loc_sys)>1)
    pos_c = peaks_pos_loc_sys(2);
    pos_d = peaks_neg_loc_sys(2);
else
    pos_c = pos_b;
    pos_d = bound_before;
end
if(~isempty(peaks_neg_loc_dias))
    pos_f = peaks_neg_loc_dias(1);
else
    pos_f = bound_after;
end

% convert positions to time points
pos_a_time = pos_a/freq;
pos_b_time = pos_b/freq;
pos_c_time = pos_c/freq;
pos_d_time = pos_d/freq;
pos_f_time = pos_f/freq;
bound_before_time = bound_before/freq;
bound_after_time = bound_after/freq;
dic_ref_time = dic_ref/freq;

%% initial values
if (isempty(x0))
    % 1. wave
    x0(1)=0.7*ppg_signal(pos_a); % amplitude (a)
    x0(2)=pos_a_time; % position (my)
    x0(3)=pos_a_time/3; % width (sigma)
    % 2. wave
    x0(4)=0.9*max([ppg_signal(pos_b) ppg_signal(pos_c) ppg_signal(pos_d)]); % amplitude (a)
    x0(5)=pos_b_time; % position (my)
    x0(6)=pos_b_time/3; % width (sigma)
    % 3. wave
    x0(7)=0.5*max([ppg_signal(pos_b) ppg_signal(pos_c) ppg_signal(pos_d)]); % amplitude (a)
    x0(8)=pos_d_time; % position (my)
    x0(9)=pos_d_time/3; % width (sigma)
    % 4. wave
    x0(10)=0.8*max(ppg_signal(bound_after:end)); % amplitude (a)
    x0(11)=pos_f_time; % position (my)
    x0(12)=min([pos_f_time (t_ppg(end)-pos_f_time)/3]); % width (sigma)
    % 5. wave
    x0(13)=0.3*max(ppg_signal(bound_after:end)); % amplitude (a)
    x0(14)=t_ppg(ppg_signal==max(ppg_signal(bound_after:end))); % position (my)
    x0(15)=t_ppg(end)-mean([t_ppg(end) pos_f_time]); % width (sigma)
    
    % i=5 in couceiro nicht konsistent definiert, da B_dias auﬂerhalb von
    % pos_f und T liegen kann; hier sind Grenzen auf bound_after und T
    % ge‰ndert worden
    
end

%% boundaries
lb=[0.5*ppg_signal(pos_a) pos_a_time 0 ...
    0.5*max([ppg_signal(pos_b) ppg_signal(pos_c) ppg_signal(pos_d)]) pos_a_time pos_a_time/3 ...
    0.2*max([ppg_signal(pos_b) ppg_signal(pos_c) ppg_signal(pos_d)]) pos_b_time pos_b_time/3 ...
    0 bound_after_time 0 ...
    0 pos_f_time 0]; % lower boundary
ub=[ppg_signal(pos_b) pos_b_time pos_b_time/3 ...
    max([ppg_signal(pos_b) ppg_signal(pos_c) ppg_signal(pos_d)]) pos_c_time pos_d_time/3 ...
    0.8*max([ppg_signal(pos_b) ppg_signal(pos_c) ppg_signal(pos_d)]) bound_before_time bound_before_time/3 ...
    max(ppg_signal(bound_after:end)) max(t_ppg) bound_after_time ...
    max(ppg_signal(bound_after:end)) max(t_ppg) bound_after_time]; % upper boundary

%% optimization
% 1. wave
gauss1=@(x) (x(1)*exp(-(t_ppg-x(2)).^2/(2*x(3)^2)));
% 2. wave
gauss2=@(x) (x(4)*exp(-(t_ppg-x(5)).^2/(2*x(6)^2)));
% 3. wave
gauss3=@(x) (x(7)*exp(-(t_ppg-x(8)).^2/(2*x(9)^2)));
% 4. wave
gauss4=@(x) (x(10)*exp(-(t_ppg-x(11)).^2/(2*x(12)^2)));
% 5. wave
gauss5=@(x) (x(13)*exp(-(t_ppg-x(14)).^2/(2*x(15)^2)));
% combining the waves
g=@(x) gauss1(x)+gauss2(x)+gauss3(x)+gauss4(x)+gauss5(x);
h=@(x) sum((ppg_signal-g(x)).^2); % Residual Sum of Squares (RSS)

% specify constraints
function [c, ceq] = constr_couceiro5kernels(x)
c = [x(1)-x(4); % 1. amplitude is less than 2. amplitude
    x(1)-x(7); % 1. amplitude is less than 3. amplitude
    x(7)-x(4); % 3. amplitude is less than 2. amplitude
    x(10)-x(4); % 4. amplitude is less than 2. amplitude
    x(13)-x(4); % 5. amplitude is less than 2. amplitude
    x(13)-x(10); % 5. amplitude is less than 4. amplitude
    x(2)-x(5); % 1. mean is less than 2. mean
    x(5)-x(8); % 2. mean is less than 3. mean
    x(8)-x(11); % 3. mean is less than 4. mean
    x(11)-x(14)]; % 4. mean is less than 5. mean
ceq = [];
end

if(~(isempty(noOpt)) && (noOpt==true)) % is noOpt is set true, skip optimization
    % initial values are given as output
    opt_val_sort = x0;
    % 1. wave
    y(1,:) = gauss1(opt_val_sort);
    % 2. wave
    y(2,:) = gauss2(opt_val_sort);
    % 3. wave
    y(3,:) = gauss3(opt_val_sort);
    % 4. wave
    y(4,:) = gauss4(opt_val_sort);
    % 5. wave
    y(5,:) = gauss5(opt_val_sort);
    % sum of basis functions
    signal_mod=sum(y); % reconstructed signal
    % produce cell array of basis functions
    y = mat2cell(y,[1 1 1 1 1],length(ppg_signal));
    
else % do optimization
    
    options = optimoptions(@fmincon,'MaxFunctionEvaluations',Inf);
    opt_values = fmincon(h,x0,[],[],[],[],lb,ub,@constr_couceiro5kernels,options);
    
    %% reconstruction
    % 1. wave
    yTmp(1,:) = gauss1(opt_values);
    % 2. wave
    yTmp(2,:) = gauss2(opt_values);
    % 3. wave
    yTmp(3,:) = gauss3(opt_values);
    % 4. wave
    yTmp(4,:) = gauss4(opt_values);
    % 5. wave
    yTmp(5,:) = gauss5(opt_values);

    % sorting the waves
    peakTmp=zeros(size(yTmp,1),1);
    for i=1:size(yTmp,1)
        try
            [~,peakTmp(i)]=find(yTmp(i,:)==max((yTmp(i,:))));
        catch
            peakTmp(i) = length(yTmp(i,:)); % take last sample as peak index if no peak exists 
        end
    end
    [~,cc]=sort(peakTmp);
    y=yTmp(cc,:);
    
    % sorting the optimized parameters
    for i=1:size(yTmp,1)
        opt_val_sort((1+(i-1)*3):(3+(i-1)*3))=...
            opt_values((1+(cc(i)-1)*3):(3+(cc(i)-1)*3));
    end
    
    % sum of basis functions
    signal_mod=sum(y); % reconstructed signal
    
    % produce cell array of basis functions
    y = mat2cell(y,[1 1 1 1 1],length(ppg_signal));
    
end

end