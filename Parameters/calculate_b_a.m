function [b_a] = calculate_b_a(PPGbeat,y,opt_params,algorithmName,freq)
% input:
% PPGbeat           ...     beat of PPG signal that is to be decomposed
% y                 ...     shapes of kernels based on optimized parameters
% opt_params        ...     optimized parameters of the kernels
% algorithmName     ...     algorithm that was used for the decomposition
% freq              ...     sampling frequency of input signal
%
% outputs:
% b_a       ...     amplitude of b wave of the second derivative of a PPG 
%                   beat over the amplitude of the a wave

second_deriv = deriv2(PPGbeat);
a = max(second_deriv); % find a
b = min(second_deriv); % find b
b_a = b/a; % calculate b over a

end