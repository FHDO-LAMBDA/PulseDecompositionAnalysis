clc
close all
clear all

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% load file containing Algorithmnames here
load('algorithmsCOMPLETE2.mat','algorithms');
% check that there is an algorithm for every entry in the list
for actualAlgorithm = 1:size(algorithms,1)
    if(~(isfile(['..\Algorithms\' 'decompose' algorithms{actualAlgorithm} '.m'])))
        errordlg('Unknown decomposition algorithm selected','Input Error','modal');
        return
    end
end
% specify source folder
sourceFolder = 'basicModels\Parameter\';
% create path for results
resultsFolder = 'results\';
caseFolder = [resultsFolder 'Cases\'];
errorFolder = [resultsFolder 'CasesError\'];
if(exist(caseFolder,'dir')~=7)
    mkdir(caseFolder)
end
if(exist(errorFolder,'dir')~=7)
    mkdir(errorFolder)
end

%% Parameters
noiseMode={'white','pink'}; % white or pink allowed
numberIterations = 25;
numberClasses = 4;
SNRvalues = [10 20 30 40 50 60 70 80 90 100];
freq = 1000; % sampling frequency
forceRevision = false; % if true snr will be calculated for all specified algorithms even if there are results from previous calculations

%% Do evaluation
for actualNoiseMode=1:2 %size(noiseMode,2) % two noise modes
    if(strcmp(noiseMode{actualNoiseMode},'white')==1 || strcmp(noiseMode{actualNoiseMode},'pink')==1)
        % if possible load previous results
        if(isfile([resultsFolder 'evaluation_' noiseMode{actualNoiseMode} '.mat']))
            load([resultsFolder 'evaluation_' noiseMode{actualNoiseMode} '.mat'],'results');
        end
        for iteration=1:numberIterations % 25 iterations
            iterationName = ['iteration' num2str(iteration)];
            for snrdb=SNRvalues % SNR in dB
                for class=1:numberClasses % four classes
                    className = ['class' num2str(class)];
                    for actualAlgorithm = 1:size(algorithms,1)
                        %% check if results already exist
                        % needs improvement for handling first run
                        %if(exist('results','var')==1 && ~forceRevision)
                        if(0)
                            if(isfield(results.(iterationName).(className),algorithms{actualAlgorithm}))
                                continue
                            end
                        end
                        %% load model signal specific for each decomposition algorithm
                        load([sourceFolder  algorithms{actualAlgorithm} '\' num2str(class) '.mat'],'signal_mod');
                        ppg = (signal_mod-min(signal_mod))/(max(signal_mod)-min(signal_mod));
                        
                        %% make noisy signal
                        noise = noise_ppg(ppg,noiseMode{actualNoiseMode},snrdb,freq);
                        
                        %% decomposition, reconstruction and calculaton of NRMSE
                        try
                            [NRMSE,~,~,~] = calculateNRMSE(noise,ppg,freq,algorithms{actualAlgorithm});
                            if(exist([caseFolder algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '\ite_' num2str(iteration) '\'],'dir')~=7)
                                mkdir([caseFolder algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '\ite_' num2str(iteration) '\'])
                            end
                            save([caseFolder algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '\ite_' num2str(iteration) '\snr_' num2str(snrdb) '_class_' num2str(class) '.mat'],'noise','ppg','freq')
                        catch
                            NRMSE = NaN;
                            if(exist([errorFolder algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '\ite_' num2str(iteration) '\'],'dir')~=7)
                                mkdir([errorFolder algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '\ite_' num2str(iteration) '\'])
                            end
                            save([errorFolder algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '\ite_' num2str(iteration) '\snr_' num2str(snrdb) '_class_' num2str(class) '.mat'],'noise','ppg','freq')
                        end
                        
                        
                        results.(iterationName).(className)(snrdb).class = class;
                        results.(iterationName).(className)(snrdb).noiseMode = noiseMode{actualNoiseMode};
                        results.(iterationName).(className)(snrdb).(algorithms{actualAlgorithm}) = NRMSE;
                        
                    end
                    results.(iterationName).(className) = orderfields(results.(iterationName).(className)); % orders fiels by name
                end
            end
        end
        save([resultsFolder 'evaluation_' noiseMode{actualNoiseMode} '.mat'],'results')
    end
end