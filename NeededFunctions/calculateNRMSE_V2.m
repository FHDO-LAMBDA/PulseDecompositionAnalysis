function [nrmse,nrmseDeriv,reconstructedSignal,kernels,opt_vals] = calculateNRMSE_V2(inputSignal,referenceSignal,freq,algorithmName)
% input:
% inputSignal           ...     signal that is to be decomposed
% referenceSignal       ...     signal that serves as a reference for NRMSE
%                               calculation
% freq                  ...     sampling frequency
% algorithmName         ...     name of the algorithm that is used for the
%                               decomposition
%
% outputs:
% nrmse                 ...     normalized root mean square error of
%                               reconstructed signal vs reference signal
% reconstructedSignal   ...     sum of kernels
% kernels               ...     optimized basis functions that result from
%                               the decomposition
% opt_vals              ...     optimized parameters of the basis functions

%% call decomposition function
decomposeFuncName = ['decompose' algorithmName];
[reconstructedSignal,kernels,opt_vals] = feval(decomposeFuncName,inputSignal,freq);

%% calculate second derivatives and NRMSE
% no derivative NRMSE
RSS = sum((referenceSignal-reconstructedSignal).^2);
STD = sum((referenceSignal-mean(referenceSignal)).^2);
nrmse=1-(sqrt(RSS)/sqrt(STD));
% derivative NRMSE
deriv_ref = deriv2(referenceSignal);
deriv_rec = deriv2(reconstructedSignal);
RSS_deriv2 = sum((deriv_ref-deriv_rec).^2);
std_deriv2 = sum((deriv_ref-mean(deriv_ref)).^2);
nrmseDeriv=1-(sqrt(RSS_deriv2)/sqrt(std_deriv2));

end