function [signal_mod,y,opt_val_sort] = decomposeGauss2(ppg_signal,freq,varargin)
% Pulse Decomposition Analysis: Gauss2 algorithm (free positioning)
%
% References:
% -
%
% This function decomposes a single beat of a ppg signal by using 2
% gaussian functions as basis functions with empirically determined initial
% values. The optimization algorithm sorts the basis functions by their 
% respective maxima location. The outputs are the 2 optimized basis 
% functions, the optimized parameters of the kernels as well as the 
% reconstructed signal which is the sum of the 2 basis functions.
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
% 'opt_val_sort':   optimized values of parameters that determine the basis
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

%% initial values
if (isempty(x0))
    % 1. wave
    x0(1)=0.8*max(ppg_signal); % amplitude (a)
    x0(2)=(2/7)*max(t_ppg); % position (my)
    x0(3)=(((2/7)*max(t_ppg))/(2*sqrt(2*log(2)))); % width (sigma)
    % 2. wave
    x0(4)=0.5*max(ppg_signal); % amplitude (a)
    x0(5)=(4/7)*max(t_ppg); % position (my)
    x0(6)=(((3/7)*max(t_ppg))/(2*sqrt(2*log(2)))); % width (sigma)
end

%% boundaries
lb=[0 0 0 0 0 0]; % lower boundary
ub=[max(ppg_signal) max(t_ppg) max(t_ppg) ...
    max(ppg_signal) max(t_ppg) max(t_ppg)]; % upper boundary

%% optimization
% 1. wave
gauss1=@(x) (x(1)*exp(-(t_ppg-x(2)).^2/(2*x(3)^2)));
% 2. wave
gauss2=@(x) (x(4)*exp(-(t_ppg-x(5)).^2/(2*x(6)^2)));
% combining the waves
g=@(x) gauss1(x)+gauss2(x);
h=@(x) sum((ppg_signal-g(x)).^2); % Residual Sum of Squares (RSS)

% specify constraints
function [c, ceq] = constr_free2kernels(x)
c = [x(2)-x(5); % 1. mean is less than 2. mean
    x(4)-x(1)]; % 1. amplitude is higher than 2. amplitude
ceq = [];
end

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
    
    options = optimoptions(@fmincon,'MaxFunctionEvaluations',Inf);
    opt_values=fmincon(h,x0,[],[],[],[],lb,ub,@constr_free2kernels,options);
    
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
