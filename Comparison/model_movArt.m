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
freq = 1000; % sampling frequency
movArt = {'mov1','mov2'}; % mov1 and mov2 allowed
numberClasses = 4;
snrdb = 100; % dummy SNR; value does not matter for this script
forceRevision = false; % if true snr will be calculated for all specified algorithms even if there are results from previous calculations

%% Do evaluation
for actualMovArt = 1:size(movArt,2)
    if(strcmp(movArt{actualMovArt},'mov1')==1 || strcmp(movArt{actualMovArt},'mov2')==1)
        % if possible load previous results
        if(isfile([resultsFolder 'evaluation_' movArt{actualMovArt} '.mat']))
            load([resultsFolder 'evaluation_' movArt{actualMovArt} '.mat'],'results');
        end
        for class=1:numberClasses
            className = ['class' num2str(class)];
            for actualAlgorithm = 1:size(algorithms,1)
                %% check if results already exist
                if(0)
                    %if(isfield(results.(className),algorithms{actualAlgorithm}) && ~forceRevision)
                    continue
                end
                %% load model signal specific for each decomposition algorithm
                load([sourceFolder  algorithms{actualAlgorithm} '\' num2str(class) '.mat'],'signal_mod');
                ppg = (signal_mod-min(signal_mod))/(max(signal_mod)-min(signal_mod));
                
                %% make noisy signal
                noise = noise_ppg(ppg,movArt{actualMovArt},snrdb,freq);
                
                %% decomposition, reconstruction and calculaton of NRMSE
                try
                    NRMSE = calculateNRMSE(noise,ppg,freq,algorithms{actualAlgorithm});
                    if(exist([caseFolder algorithms{actualAlgorithm} '\'],'dir')~=7)
                        mkdir([caseFolder algorithms{actualAlgorithm} '\'])
                    end
                    save([caseFolder algorithms{actualAlgorithm} '\noise_' movArt{actualMovArt} 'class_' num2str(class) '.mat'],'noise','ppg','freq')
                catch
                    NRMSE = NaN;
                    if(exist([errorFolder algorithms{actualAlgorithm} '\'],'dir')~=7)
                        mkdir([errorFolder algorithms{actualAlgorithm} '\'])
                    end
                    save([errorFolder algorithms{actualAlgorithm} '\noise_' movArt{actualMovArt} 'class_' num2str(class) '.mat'],'noise','ppg','freq')
                end
                
                results.(className).(algorithms{actualAlgorithm}) = NRMSE;
                
            end
            results.(className) = orderfields(results.(className)); % orders fiels by name
        end
        save([resultsFolder 'evaluation_' movArt{actualMovArt} '.mat'],'results')
    end
end