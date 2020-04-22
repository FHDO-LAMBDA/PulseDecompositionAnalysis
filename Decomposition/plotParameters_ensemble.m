clear all
close all
clc
%
% takes only successful decompositions into account
%

%% Paths and initialization
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% add path to parameter functions and cell with parameter names
addpath('..\Parameters');
% load file containing parameters to be plotted here
load('parameters.mat');
% load file containing algorithms here
load('algorithmsCOMPLETE.mat','algorithms');
% check that there is an algorithm for every entry in the list
for actualAlgorithm = 1:size(algorithms,1)
    if(~(isfile(['..\Algorithms\' 'decompose' algorithms{actualAlgorithm} '.m'])))
        errordlg('Unknown decomposition algorithm selected','Input Error','modal');
        return
    end
end
% specify plot, source and results folder
sourceFolder='realData\';
dataFolder='decompositionParameters_ensemble\';
plotFolder='plots_ensemble\';
% load data information
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']);
patients=physiologicalMeasuresTable.SubjectID; % loads list with patients
dataset = cell(size(algorithms,1),size(epochs,1),size(parameters,1));

%% Format data
for actualPatientNumber = 1:size(patients,1)
    fileID = patients{actualPatientNumber};
    for currentInterval = 1:size(epochs,1)
        if(exist([dataFolder fileID '\' epochs{currentInterval}],'dir') ~= 7)
            continue
        end
        for actualAlgorithm = 1:size(algorithms,1)
            for actualParameter = 1:size(parameters,1)
                try
                    load([dataFolder fileID '\' epochs{currentInterval} '\' algorithms{actualAlgorithm} '.mat'],parameters{actualParameter,1});
                    dataset{actualAlgorithm,currentInterval,actualParameter}(end+1) = eval(parameters{actualParameter,1});
                    clear(parameters{actualParameter,1});
                catch
                    dataset{actualAlgorithm,currentInterval,actualParameter}(end+1) = NaN;
                end
            end
        end
    end
end

matrix = zeros(size(algorithms,1),size(epochs,1),size(parameters,1),size(patients,1));
matrix(matrix==0) = NaN; % initialize matrix

% fill matrix (for boxplot comparing epochs for each algorithm)
for actualAlgorithm = 1:size(algorithms,1)
    for currentInterval = 1:size(epochs,1)
        for actualParameter = 1:size(parameters,1)
            try
                matrix(actualAlgorithm,currentInterval,actualParameter,1:length(dataset{actualAlgorithm,currentInterval,actualParameter}(~isnan(dataset{actualAlgorithm,currentInterval,actualParameter})))) = dataset{actualAlgorithm,currentInterval,actualParameter}(~isnan(dataset{actualAlgorithm,currentInterval,actualParameter}));
            catch
                matrix(actualAlgorithm,currentInterval,actualParameter,1:length(dataset{actualAlgorithm,currentInterval,actualParameter}(~isnan(dataset{actualAlgorithm,currentInterval,actualParameter})))) = NaN;
            end
        end
    end
end

% reshape matrix (for boxplot comparing the algorithms)
matrix_ges = zeros(size(matrix,2)*size(matrix,4),size(matrix,1),size(matrix,3));
for actualAlgorithm = 1:size(algorithms,1)
    for currentInterval = 1:size(epochs,1)
        for actualParameter = 1:size(parameters,1)
            matrix_ges(1+((currentInterval-1)*size(matrix,4)):currentInterval*size(matrix,4),actualAlgorithm,actualParameter) = matrix(actualAlgorithm,currentInterval,actualParameter,:);
        end
    end
end

%% Do actual plotting
% comparison of allgorithms as a whole
for actualParameter = 1:size(parameters,1)
    myAxis = [-inf inf parameters{actualParameter,2} parameters{actualParameter,3}];
    figure('Name',parameters{actualParameter,1});
    boxplot(matrix_ges(:,:,actualParameter),'Labels',algorithms')
    ylabel(parameters{actualParameter,1},'Interpreter','none')
    xlabel('algorithms')
    xtickangle(45)
    axis(myAxis)
    currentFigure=gcf;
    if(exist([plotFolder parameters{actualParameter,1} '\'],'dir')~=7)
        mkdir([plotFolder parameters{actualParameter,1} '\'])
    end
    savefig([plotFolder parameters{actualParameter,1} '\' parameters{actualParameter,1} '.fig']);
    saveas(currentFigure,[plotFolder parameters{actualParameter,1} '\' parameters{actualParameter,1} '.png']);
    print(currentFigure,[plotFolder parameters{actualParameter,1} '\' parameters{actualParameter,1}],'-dpdf');
    matlab2tikz([plotFolder parameters{actualParameter,1} '\' parameters{actualParameter,1} '.tex']);
end

% comparison of epochs for each algorithm
for actualAlgorithm = 1:size(algorithms,1)
    for actualParameter = 1:size(parameters,1)
        myAxis = [-inf inf parameters{actualParameter,2} parameters{actualParameter,3}];
        currentMatrix = squeeze(matrix(actualAlgorithm,:,actualParameter,:));
        figure('Name',algorithms{actualAlgorithm});
        boxplot(currentMatrix','Labels',epochs')
        ylabel(parameters{actualParameter,1},'Interpreter','none')
        xlabel('epochs')
        xtickangle(45)
        axis(myAxis)
        currentFigure=gcf;
        if(exist([plotFolder parameters{actualParameter,1} '\'],'dir')~=7)
            mkdir([plotFolder parameters{actualParameter,1} '\'])
        end
        savefig([plotFolder parameters{actualParameter,1} '\' algorithms{actualAlgorithm,1} '.fig']);
        saveas(currentFigure,[plotFolder parameters{actualParameter,1} '\' algorithms{actualAlgorithm,1} '.png']);
        print(currentFigure,[plotFolder parameters{actualParameter,1} '\' algorithms{actualAlgorithm,1}],'-dpdf');
        matlab2tikz([plotFolder parameters{actualParameter,1} '\' algorithms{actualAlgorithm,1} '.tex']);
    end
end
