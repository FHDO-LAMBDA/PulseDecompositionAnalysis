clear all
close all
clc

% takes only successful decompositions into account
% muss noch angepasst werden an mehrere Datensätze


data ='CPT';
sourceFolder=['Data\' data '\realData\'];
resultsFolder=['Data\' data '\results_ensembleReference\'];
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']); 
load('algorithmsCOMPLETE2.mat'); 
patients=physiologicalMeasuresTable.SubjectID;
dataset = cell(size(algorithms,1),size(epochs,1));
datasetD2 = cell(size(algorithms,1),size(epochs,1));

%% format data
for actualPatientNumber = 1:size(patients,1)
    fileID = patients{actualPatientNumber};
    for currentInterval = 1:3 % size(epochs,1) % 1:3 ist Baseline
        if(exist([resultsFolder fileID '\' epochs{currentInterval}]) ~= 7)
            continue
        end
        for actualAlgorithm = 1:size(algorithms,1)
            try
                load([resultsFolder fileID '\' epochs{currentInterval} '\' algorithms{actualAlgorithm} '.mat'],'nrmse','nrmseD2');
                dataset{actualAlgorithm,currentInterval}(end+1) = nrmse;
                datasetD2{actualAlgorithm,currentInterval}(end+1) = nrmseD2;
                clear('nrmse','nrmseD2');
            catch
                dataset{actualAlgorithm,currentInterval}(end+1) = NaN;
                datasetD2{actualAlgorithm,currentInterval}(end+1) = NaN;
            end
        end
    end
end

matrix = zeros(size(algorithms,1),size(epochs,1),size(patients,1));
matrixD2 = zeros(size(algorithms,1),size(epochs,1),size(patients,1));
matrix(matrix==0) = NaN; % initialize matrix
matrixD2(matrixD2==0) = NaN; % initialize matrix

% fill matrix (for boxplot comparing epochs for each algorithm)
for actualAlgorithm = 1:size(algorithms,1)
    for currentInterval = 1:size(epochs,1)

            try
                matrix(actualAlgorithm,currentInterval,1:length(dataset{actualAlgorithm,currentInterval}(~isnan(dataset{actualAlgorithm,currentInterval})))) = dataset{actualAlgorithm,currentInterval}(~isnan(dataset{actualAlgorithm,currentInterval}));
            catch
                matrix(actualAlgorithm,currentInterval,1:length(dataset{actualAlgorithm,currentInterval}(~isnan(dataset{actualAlgorithm,currentInterval})))) = NaN;
            end
            
            try
                matrixD2(actualAlgorithm,currentInterval,1:length(datasetD2{actualAlgorithm,currentInterval}(~isnan(datasetD2{actualAlgorithm,currentInterval})))) = datasetD2{actualAlgorithm,currentInterval}(~isnan(datasetD2{actualAlgorithm,currentInterval}));
            catch
                matrixD2(actualAlgorithm,currentInterval,1:length(datasetD2{actualAlgorithm,currentInterval}(~isnan(datasetD2{actualAlgorithm,currentInterval})))) = NaN;
            end

    end
end

% reshape matrix (for boxplot comparing the algorithms)
matrix_ges = zeros(size(matrix,2)*size(matrix,3),size(matrix,1));
matrix_gesD2 = zeros(size(matrixD2,2)*size(matrixD2,3),size(matrixD2,1));
for actualAlgorithm = 1:size(algorithms,1)
    for currentInterval = 1:size(epochs,1)
        matrix_ges(1+((currentInterval-1)*size(matrix,3)):currentInterval*size(matrix,3),actualAlgorithm) = matrix(actualAlgorithm,currentInterval,:);
        matrix_gesD2(1+((currentInterval-1)*size(matrixD2,3)):currentInterval*size(matrixD2,3),actualAlgorithm) = matrixD2(actualAlgorithm,currentInterval,:);
    end
end

%% make pairwise comparison of algorithms
% determine pValues
pMatrix = zeros(size(algorithms,1));
pTable = table('Size',[numel(algorithms) numel(algorithms)],...
            'VariableTypes',repmat({'double'},[1 numel(algorithms)]),...
            'VariableNames',algorithms,...
            'RowNames',algorithms);
for i = 1:numel(algorithms)
    for j = 1:numel(algorithms)
        [~,pValue] = ttest2(matrix_gesD2(:,i),matrix_gesD2(:,j));
        pMatrix(i,j) = pValue;
        pTable{i,j} = pValue;
    end
end

% correct pValues
pMatrixCorrected = zeros(size(algorithms,1));
pTableCorrected = table('Size',[numel(algorithms) numel(algorithms)],...
            'VariableTypes',repmat({'double'},[1 numel(algorithms)]),...
            'VariableNames',algorithms,...
            'RowNames',algorithms);
alpha = 0.05;
pPostHoc = zeros(nchoosek(numel(algorithms),2),1);
counter = 1;
for i=1:numel(algorithms)
    for j=1:numel(algorithms)
        if(i < j)
            pPostHoc(counter) = pTable{i,j};
            counter = counter + 1;
        end
    end
end
[pPostHocCorrected,~] = bonf_holm(pPostHoc,alpha);
for currentPValue = 1:numel(pPostHoc)
%     if(pPostHocCorrected(currentPValue)>1)
%        pPostHocCorrected(currentPValue) = 1; 
%     end
    pMatrixCorrected(pMatrix==pPostHoc(currentPValue)) = pPostHocCorrected(currentPValue);
end
for i=1:numel(algorithms)
    for j=1:numel(algorithms)
        if(i==j)
            pTableCorrected{i,j} = 1;
        else
            pTableCorrected{i,j} = pMatrixCorrected(i,j);
        end
    end
end

% save tables (corrected, uncorrected)
save([resultsFolder 'pTable.mat'],'pTable','pTableCorrected')

% save a textfile that can be inserted in latex and contains values for
% table
fileID = fopen([resultsFolder 'tableLatex.txt'],'w');
for i = 1:numel(algorithms)
    fprintf(fileID,'\\textbf{%d} & ',i);
    for j = 1:numel(algorithms)-1
        if(pTableCorrected{i,j}<0.05)
            fprintf(fileID,'\\textcolor{blue}{%.2f} & ',pTableCorrected{i,j});
        else
            fprintf(fileID,'n.s. & ');
        end
    end
    if(pTableCorrected{i,numel(algorithms)}<0.05)
        fprintf(fileID,'\\textcolor{blue}{%.2f}',pTableCorrected{i,numel(algorithms)});
    else
        fprintf(fileID,'n.s.');
    end
    fprintf(fileID,'\\\\ \n');
    if(~(i==numel(algorithms)))
        fprintf(fileID,'\\mr \n');
    end
end
fclose(fileID);

%% make figures
myAxis = [-inf inf 0 1];

% comparison of allgorithms as a whole
figure('Name','nrmse');
boxplot(matrix_ges,'Labels',algorithms','Colors','k','Widths',0.3);
h = findobj(gca,'tag','Outliers');
delete(h) % delete outliers from plot
ylabel('NRMSE','Interpreter','none')
xtickangle(45)
axis(myAxis)
set(gcf, 'Units', 'centimeters');
set(gcf, 'PaperUnits', 'centimeters');
currentPos=get(gcf,'Position');
set(gcf, 'Position', [currentPos(1) currentPos(1) 18 8]);
set(findobj(gcf,'type','axes'),...
    'FontSize', 9,...
    'FontName', 'Times New Roman',...
    'FontWeight','normal',...
    'TitleFontWeight','normal');
savefig([resultsFolder 'NRMSE.fig'])
matlab2tikz([resultsFolder 'NRMSE.tex'])

% comparison of allgorithms as a whole
figure('Name','nrmse_abl');
boxplot(matrix_gesD2,'Labels',algorithms','Colors','k','Widths',0.3);
h = findobj(gca,'tag','Outliers');
delete(h) % delete outliers from plot
ylabel('NRMSE','Interpreter','none')
xtickangle(45)
axis(myAxis)
set(gcf, 'Units', 'centimeters');
set(gcf, 'PaperUnits', 'centimeters');
currentPos=get(gcf,'Position');
set(gcf, 'Position', [currentPos(1) currentPos(1) 18 8]);
set(findobj(gcf,'type','axes'),...
    'FontSize', 9,...
    'FontName', 'Times New Roman',...
    'FontWeight','normal',...
    'TitleFontWeight','normal');
savefig([resultsFolder 'NRMSEabl.fig'])
matlab2tikz([resultsFolder 'NRMSEabl.tex'])
