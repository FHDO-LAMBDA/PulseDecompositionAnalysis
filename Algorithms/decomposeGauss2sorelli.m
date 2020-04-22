function [signal_mod,y,opt_val_sort] = decomposeGauss2sorelli(ppg_signal,freq,varargin)
% Pulse Decomposition Analysis: Gauss2 algorithm (Sorelli restrictions)
%
% References:
% Michele Sorelli, Antonia Perrella und Leonardo Bocchi. „Detecting 
% vascular age using the analysis of peripheral pulse“. In: IEEE 
% Transactions on Biomedical Engineering PP.c (2018), S. 1.
%
% NOTE:         This version works as close to the description in Sorellis 
%               paper as possible.
% Differences:  - forward difference instead of central difference is used
%                 for derivation
%
% This function decomposes a single beat of a ppg signal by using 2 gaussian
% functions as basis functions with the restrictions of Sorelli et al. The 
% optimization algorithm sorts the basis functions by their respective 
% maxima location. The outputs are the 2 optimized basis functions, the 
% optimized parameters of the kernels as well as the reconstructed signal 
% which is the sum of the 2 basis functions.
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
% 'y':              cell array of size 2x1 which contains the 2 basis
%                   functions (y1, y2) each with a size of 
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

%% Sorellis algorithm - finding initial values
if (isempty(x0))
    % 1. search for systolic reference (highest maximum)
    [peaks,~] = findpeaks(ppg_signal);
    sys_ref = find(ppg_signal==max(peaks));
    
    % 2. find end-diastolic perfusion through as absolute minimum (but not
    % first point in signal)
    leave_out = 5;
    dia_ref = find(ppg_signal(leave_out+1:end)==min(ppg_signal(leave_out+1:end))); % leave out few first points
    dia_ref = dia_ref+leave_out; % add samples to be in line with real signal input
    dia_ref = t_ppg(dia_ref); % turn diastolic reference to time point instead of sample
    
    % 3. refinement step for systolic detection: search for all maxima
    % preceding the original solution, whose prominence with respect to the
    % corresponding valley does not fall below 80% of the original peak
    % amplitude
    [~,maxima_loc] = findpeaks(ppg_signal,'MinPeakProminence',0.8*ppg_signal(sys_ref)); % find all peaks with required prominence
    maxima_loc(maxima_loc>=sys_ref) = []; % keep only maxima before systolic reference
    if ~(isempty(maxima_loc))
        sys_ref = maxima_loc(1); % make first found maximum to new systolic reference
        % this step is not explicitly stated in sorellis Paper
    end
    sys_ref = t_ppg(sys_ref); % turn systolic reference to time point instead of sample
    
    % 4. refinement step for diastolic detection: not needed in this function
    % because input signals should be beats from minimum to minimum
    
    % 5. calculte derivative (by 3point differentiator) (zentrierte differenz)
    deriv_1 = deriv1(ppg_signal); % not using central difference
    
    % 6. filter derivative with 7point moving average
    windowWidth = 7;
    b = (1/windowWidth)*ones(1,windowWidth);
    a = 1;
    deriv_1 = filtfilt(b,a,deriv_1);
    
    % 7. search for earliest neg to pos zero crossing of first derivative after
    % systolic reference. If none are detected, end-diastolic reference is
    % selected
    sign_deriv_1 = sign(deriv_1); % Array mit Vorzeichen (-1 = negativ, 1= positiv)
    zero_cross = t_ppg([diff(sign_deriv_1) 0] ~= 0); % find time stamps of sign change
    if(~(isempty(zero_cross)))
        zero_cross(zero_cross<=sys_ref) = [];
        if(~(isempty(zero_cross)))
            cross_ref = zero_cross(1);
        else
            cross_ref = dia_ref;
        end
    else
        cross_ref = dia_ref;
    end
    
    % 8. time span between this point and the systolic peak is analyzed for the
    % presence of local p'(t) maxima exceeding the average pulse slope (average
    % of first derivative in this interval) in the same interval: if detected,
    % the earliest of them is adopted as the incisura reference, otherwise the
    % original p'(t) zero crossing is chosen
    av_slope = mean(deriv_1(find(t_ppg==sys_ref):find(t_ppg==cross_ref)));
    try
        [~,inc_ref] = findpeaks(deriv_1(find(t_ppg==sys_ref):find(t_ppg==cross_ref)),'MinPeakHeight',av_slope);
    catch
        inc_ref = [];
    end
    if(~(isempty(inc_ref)))
        inc_ref = inc_ref+(sys_ref*freq); % am ende anzahl samples von sys_ref wieder raufrechnen
        inc_ref = inc_ref/freq; % turn incisura reference to time point instead of sample
    else
        if(~(cross_ref==dia_ref))
            inc_ref = cross_ref;
        else
            inc_ref = sys_ref;
        end
    end
    
    %% calculate systolic and diastolic time span
    inc_ref = inc_ref(1); % make sure there is only one inc_ref
    delta_t_sys = inc_ref; % delta_t_sys is time from beginning to incisura
    delta_t_dia = t_ppg(end)-inc_ref; % delta_t_sys is time from incisura to the end
    
    %% initial values
    % 1. wave
    x0(1)=0.8*max(ppg_signal); % amplitude (a)
    x0(2)=t_ppg(1)+0.5*delta_t_sys; % position (my)
    x0(3)=delta_t_sys/(2*sqrt(2*log(2))); % width (sigma)
    % 2. wave
    x0(4)=0.4*max(ppg_signal); % amplitude (a)
    x0(5)=inc_ref+0.33*delta_t_dia; % position (my)
    x0(6)=0.75*delta_t_dia/(2*sqrt(2*log(2))); % width (sigma)
end

%% boundaries
lb=[0.5*max(ppg_signal) t_ppg(1) 0.5*delta_t_sys/(2*sqrt(2*log(2)))... 
    0 inc_ref 0.3*delta_t_dia/(2*sqrt(2*log(2)))]; % lower boundary
ub=[max(ppg_signal) inc_ref 1.5*delta_t_sys/(2*sqrt(2*log(2)))...
    0.6*max(ppg_signal) max(t_ppg) 1*delta_t_dia/(2*sqrt(2*log(2)))]; % upper boundary

%% optimization
% 1. wave
gauss1=@(x) (x(1)*exp(-((t_ppg-x(2)).^2)/(2*(x(3)^2))));
% 2. wave
gauss2=@(x) (x(4)*exp(-((t_ppg-x(5)).^2)/(2*(x(6)^2))));
% combining the waves
g=@(x) gauss1(x)+gauss2(x);
h=@(x) sum((ppg_signal-g(x)).^2); % Residual Sum of Squares (RSS)

if(~(isempty(noOpt)) && (noOpt==true)) % is noOpt is set true, skip optimization
    % initial values are given as output
    opt_val_sort = x0;
    % 1. wave
    y(1,:) = gauss1(opt_val_sort);
    % 2. wave
    y(2,:) = gauss2(opt_val_sort);
    % sum of basis functions
    signal_mod=sum(y); % reconstructed signal
    % produce cell array of basis functions
    y = mat2cell(y,[1 1],length(ppg_signal));
    
else % do optimization
    
    options = optimoptions(@fmincon,'MaxFunEvals',Inf);
    opt_values =fmincon(h,x0,[],[],[],[],lb,ub,[],options); % two-kernel sorelli algorithm does not need position constraints
    
    %% reconstruction
    % 1. wave
    yTmp(1,:) = gauss1(opt_values);
    % 2. wave
    yTmp(2,:) = gauss2(opt_values);

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
    y = mat2cell(y,[1 1],length(ppg_signal));
    
end

end