clear all
close all
clc

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% load file containing algorithms here
load('algorithmsSORELLI.mat','algorithms');
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
resultsFolder=['Data\' dataset '\results_beatwise\'];
% load data information
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']);
patients=physiologicalMeasuresTable.SubjectID; % loads list with patients

%% Parameters
samplingFreq = 100;
cameraResolution = 4096;

%% Do decomposition
for actualPatientNumber=1:size(patients,1)
    for currentInterval=1:size(epochs,1)
        if(exist([sourceFolder patients{actualPatientNumber} '\' epochs{currentInterval} '.mat'],'file') ~= 2)
            continue
        end
        %% load signal to be decomposed
        load([sourceFolder patients{actualPatientNumber} '\' epochs{currentInterval} '.mat'],...
            'beatIndices','beatIndicesEnsembleBeat','ppgGreen');
        if((~any(isnan(beatIndices))) && (~any(isnan(beatIndicesEnsembleBeat))))
            ppgGreen = cameraResolution - ppgGreen;
            singleBeats = createSingleBeats(ppgGreen,samplingFreq,beatIndices,beatIndicesEnsembleBeat);
        else
            singleBeats = NaN;
        end
        for actualAlgorithm = 1:size(algorithms,1)
            %% decomposition, reconstruction and calculation of NRMSE
            decompositionResults = struct;
            for beatNumber = 1:size(beatIndicesEnsembleBeat,1)
                if(iscell(singleBeats))
                    decompositionResults(beatNumber).singleBeats = singleBeats{beatNumber};
                elseif(isnan(singleBeats))
                    decompositionResults(beatNumber).singleBeats = NaN;
                else
                    errordlg('Unknown Error');
                    return
                end
                try
                    [nrmse,signal_mod,y,opt_params] = calculateNRMSE(singleBeats{beatNumber},singleBeats{beatNumber},samplingFreq,algorithms{actualAlgorithm});
                    decompositionResults(beatNumber).nrmse = nrmse;
                    decompositionResults(beatNumber).signal_mod = signal_mod;
                    decompositionResults(beatNumber).y = y;
                    decompositionResults(beatNumber).opt_params = opt_params;
                catch
                    decompositionResults(beatNumber).nrmse = NaN;
                    decompositionResults(beatNumber).signal_mod = NaN;
                    decompositionResults(beatNumber).y = cell(3,1);
                    decompositionResults(beatNumber).opt_params = NaN;
                end
            end
            if(exist([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\'],'dir')~=7)
                mkdir([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\'])
            end
            save([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\' algorithms{actualAlgorithm} '.mat'],...
                'decompositionResults', ...
                'beatIndices','beatIndicesEnsembleBeat', 'ppgGreen', 'samplingFreq');
            clear decompositionResults
        end
        clear ppgGreen beatIndices beatIndicesEnsembleBeat
    end
end