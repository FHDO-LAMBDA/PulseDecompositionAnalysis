clear all
close all
clc

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% add path to parameter functions and cell with parameter names
addpath('..\Parameters');
% load file containing algorithms here
load('algorithmsSORELLI.mat','algorithms');
% check that there is an algorithm for every entry in the list
for actualAlgorithm = 1:size(algorithms,1)
    if(~(isfile(['..\Algorithms\' 'decompose' algorithms{actualAlgorithm} '.m'])))
        errordlg('Unknown decomposition algorithm selected','Input Error','modal');
        return
    end
end
% specify source, data and results folder
sourceFolder='realData\';
resultsFolder='decompositionParameters_beatwise\';
dataFolder='results_beatwise\';
% load data information
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']);
patients=physiologicalMeasuresTable.SubjectID; % loads list with patients

%% Extract parameters
for actualPatientNumber=1:size(patients,1)
    for actualAlgorithm = 1:size(algorithms,1)
        for currentInterval=1:size(epochs,1)
            currentFile = [dataFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\' algorithms{actualAlgorithm}];
            % try loading current file
            try
                load(currentFile)
            catch
                b_a = NaN;
                T_sys_dia = NaN;
                continue
            end
            % initialize results struct
            decompositionParameters = struct;
            for beatNumber = 1:size(beatIndicesEnsembleBeat)
                %% calculate parameters
                % second derivative
                try
                    decompositionParameters(beatNumber).b_a = ...
                        calculate_b_a(decompositionResults(beatNumber).signal_mod,decompositionResults(beatNumber).opt_params,algorithms{actualAlgorithm});
                catch
                    decompositionParameters(beatNumber).b_a = NaN;
                end
                % kernel based
                try
                    decompositionParameters(beatNumber).T_sys_dia = ...
                        calculate_T_sys_dia(decompositionResults(beatNumber).signal_mod,decompositionResults(beatNumber).opt_params,algorithms{actualAlgorithm});
                catch
                    decompositionParameters(beatNumber).T_sys_dia = NaN;
                end
            end
            %% calculate mean of parameters
            % second derivative
            b_a = mean([decompositionParameters.b_a],'omitnan');
            % kernel based
            T_sys_dia = mean([decompositionParameters.T_sys_dia],'omitnan');
            
            %% save parameters
            if(exist([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\'],'dir')~=7)
                mkdir([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\'])
            end
            save([resultsFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\' algorithms{actualAlgorithm} '.mat'],...
                'decompositionParameters','b_a','T_sys_dia');
            clear decompositionParameters b_a T_sys_dia
            
        end
    end
end
