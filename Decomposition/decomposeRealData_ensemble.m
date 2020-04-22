clear all
close all
clc

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% load file containing algorithms here
load('algorithmsBEST.mat','algorithms');
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
resultsFolder=['Data\' dataset '\results_ensemble\'];
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
        %% load signal to be decomposed
        load([sourceFolder patients{actualPatientNumber} '\' epochs{currentInterval} '.mat'],...
            'beatStartIndex','beatStopIndex','meanSeg','beatIndicesEnsembleBeat');
        for actualAlgorithm = 1:size(algorithms,1)
            %% decomposition, reconstruction and calculation of NRMSE
            try
                [nrmse,signal_mod,y,opt_params] = calculateNRMSE(meanSeg(beatStartIndex:beatStopIndex),meanSeg(beatStartIndex:beatStopIndex),samplingFreq,algorithms{actualAlgorithm});
            catch
                nrmse = NaN;
                signal_mod = NaN;
                y = cell(3,1);
                opt_params = NaN;
            end
            if(exist([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\'],'dir')~=7)
                mkdir([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\'])
            end
            save([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\' algorithms{actualAlgorithm} '.mat'],...
                'nrmse', 'signal_mod','y', 'opt_params', ...
                'beatStartIndex', 'beatStopIndex', 'meanSeg', ...
                'beatIndicesEnsembleBeat', 'samplingFreq');
            clear signal_mod y opt_params nrmse
        end
        clear meanSeg beatStartIndex beatStopIndex beatIndicesEnsembleBeat
    end
end