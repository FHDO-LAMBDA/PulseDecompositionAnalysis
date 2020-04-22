clear all
close all
clc

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% load file containing algorithms here
load('algorithmsCOMPLETE2.mat','algorithms');
% check that there is an algorithm for every entry in the list
for actualAlgorithm = 1:size(algorithms,1)
    if(~(isfile(['..\Algorithms\' 'decompose' algorithms{actualAlgorithm} '.m'])))
        errordlg('Unknown decomposition algorithm selected','Input Error','modal');
        return
    end
end
% specify source and results folder
dataset ='CPT';
sourceFolder=['Data\' dataset '\realData\'];
resultsFolder=['Data\' dataset '\results_ensembleReference\'];
% load data information
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']);
patients=physiologicalMeasuresTable.SubjectID; % loads list with patients

%% Parameters
samplingFreq = 100;

%% Do decomposition
for actualPatientNumber=1:size(patients,1)
    for currentInterval=1:size(epochs,1)
        if(exist([sourceFolder patients{actualPatientNumber} '\' epochs{currentInterval} '.mat'],'file') ~= 2)
            continue
        end
        %% load indices
        load([sourceFolder patients{actualPatientNumber} '\' epochs{currentInterval} '.mat'],...
            'beatStartIndexReference','beatStopIndexReference');
        for actualAlgorithm = 1:size(algorithms,1)
            %% decomposition, reconstruction and calculation of NRMSE
            try
                [nrmse,nrmseD2,signal_mod,y,opt_params] = calculateNRMSE_V2(physiologicalMeasuresTable.([epochs{currentInterval} '_meanSeg']){actualPatientNumber,1}(beatStartIndexReference:beatStopIndexReference),...
                    physiologicalMeasuresTable.([epochs{currentInterval} '_meanSeg']){actualPatientNumber,1}(beatStartIndexReference:beatStopIndexReference),samplingFreq,algorithms{actualAlgorithm});
            catch
                nrmse = NaN;
                nrmseD2 = NaN;
                signal_mod = NaN;
                y = cell(3,1);
                opt_params = NaN;
            end
            if(exist([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\'],'dir')~=7)
                mkdir([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\'])
            end
            save([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\' algorithms{actualAlgorithm} '.mat'],...
                'nrmse', 'nrmseD2', 'signal_mod','y', 'opt_params', ...
                'beatStartIndexReference', 'beatStopIndexReference', 'samplingFreq');
            clear signal_mod y opt_params nrmse nrmseD2
        end
        clear beatStartIndexReference beatStopIndexReference
    end
end