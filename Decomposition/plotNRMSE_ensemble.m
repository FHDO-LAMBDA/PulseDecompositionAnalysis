clear all
close all
clc

% takes only successful decompositions into account


sourceFolder='realData\';
resultsFolder='results_ensemble\';
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']); 
load('algorithmsSORELLI.mat'); 
patients=physiologicalMeasuresTable.SubjectID;
dataset = cell(size(algorithms,1),size(epochs,1));


for actualPatientNumber = 1:size(patients,1)
    fileID = patients{actualPatientNumber};
    for currentInterval = 1:size(epochs,1)
        if(exist([resultsFolder fileID '\' epochs{currentInterval}]) ~= 7)
            continue
        end
        for actualAlgorithm = 1:size(algorithms,1)
            try
                load([resultsFolder fileID '\' epochs{currentInterval} '\' algorithms{actualAlgorithm} '.mat'],'nrmse');
                dataset{actualAlgorithm,currentInterval}(end+1) = nrmse;
                clear('nrmse');
            catch
                dataset{actualAlgorithm,currentInterval}(end+1) = NaN;
            end
        end
    end
end

matrix = zeros(size(algorithms,1),size(epochs,1),size(patients,1));
matrix(matrix==0) = NaN; % initialize matrix

% fill matrix (for boxplot comparing epochs for each algorithm)
for actualAlgorithm = 1:size(algorithms,1)
    for currentInterval = 1:size(epochs,1)

            try
                matrix(actualAlgorithm,currentInterval,1:length(dataset{actualAlgorithm,currentInterval}(~isnan(dataset{actualAlgorithm,currentInterval})))) = dataset{actualAlgorithm,currentInterval}(~isnan(dataset{actualAlgorithm,currentInterval}));
            catch
                matrix(actualAlgorithm,currentInterval,1:length(dataset{actualAlgorithm,currentInterval}(~isnan(dataset{actualAlgorithm,currentInterval})))) = NaN;
            end

    end
end

% reshape matrix (for boxplot comparing the algorithms)
matrix_ges = zeros(size(matrix,2)*size(matrix,3),size(matrix,1));
for actualAlgorithm = 1:size(algorithms,1)
    for currentInterval = 1:size(epochs,1)
        matrix_ges(1+((currentInterval-1)*size(matrix,3)):currentInterval*size(matrix,3),actualAlgorithm) = matrix(actualAlgorithm,currentInterval,:);
    end
end
% reshape does not place elements right
%matrix_ges = reshape(matrix,[size(matrix,2)*size(matrix,4),size(matrix,1),size(matrix,3)]);

% make figures
myAxis = [-inf inf -1 1];

% comparison of allgorithms as a whole
figure('Name','nrmse_abl');
boxplot(matrix_ges,'Labels',algorithms')
ylabel('nrmse_abl','Interpreter','none')
xlabel('algorithms')
xtickangle(45)
axis(myAxis)


% comparison of epochs for each algorithm
for actualAlgorithm = 1:size(algorithms,1)
    currentMatrix = squeeze(matrix(actualAlgorithm,:,:));
    figure('Name',algorithms{actualAlgorithm});
    boxplot(currentMatrix','Labels',epochs')
    ylabel('nrmse_abl','Interpreter','none')
    xlabel('epochs')
    xtickangle(45)
    axis(myAxis)
end
