function [RI_peaks] = calculate_RI_peaks(PPGbeat,y,opt_params,algorithmName,freq)
% input:
% PPGbeat           ...     beat of PPG signal that is to be decomposed
% y                 ...     shapes of kernels based on optimized parameters
% opt_params        ...     optimized parameters of the kernels
% algorithmName     ...     algorithm that was used for the decomposition
% freq              ...     sampling frequency of input signal
%
% outputs:
% RI_peaks          ...     amplitude of systolic component over amplitude
%                           of diastolic component

%% exceptions
% GammaGauss4
if(strcmp(algorithmName,'GammaGauss4'))
    RI_peaks = NaN;
    return
end

% opt_params is NaN
if(isnan(opt_params))
    RI_peaks = NaN;
    return
end

%% calculate ratio of amplitudes of systolic and diastolic component
numKernels = length(opt_params)/3; % get number of kernels
switch numKernels
    case 0
        errordlg('Number of kernels is 0. No calculation of area ratio possible','Input Error','modal');
        return
    case 1
        errordlg('Number of kernels is 1. No calculation of area ratio possible','Input Error','modal');
        return
    case 2
        curve_sys = y{1};
        P_sys = max(curve_sys);
        curve_dia = y{2};
        P_dia = max(curve_dia);
    case 3
        curve_sys = y{1};
        P_sys = max(curve_sys);
        curve_dia = sum([y{2};y{3}]);
        P_dia = max(curve_dia);
    case 4
        curve_sys = y{1};
        P_sys = max(curve_sys);
        curve_dia = sum([y{2};y{3};y{4}]);
        P_dia = max(curve_dia);
    case 5
        curve_sys = sum([y{1};y{2}]);
        P_sys = max(curve_sys);
        curve_dia = sum([y{3};y{4};y{5}]);
        P_dia = max(curve_dia);
    otherwise
        error('Number of kernels exceeds 5.');
end
RI_peaks = P_dia/P_sys;

end