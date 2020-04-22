%% TODO: möglichkeit einführen, nur bestimmte Paarungen in post hoc zu testen

% in der physiological measures table habe ich sogar das referenz PPG
% vllt sollte ich die Zerlegung auch in einer Tabelle abspeichern?

% abfragen, ob überhaupt für alle algorithmen, die ich spezifiziere eine
% Tabelle da ist und Werte für alle Parameter?

% auswahl, ob Korrektur überhaupt nötig bei omnibus test (oder ist die
% immer nötig bei mehr als einem Parameter?)

% check ob epochsofinteresr etc teilmenge von epochs sind. bei Parametern
% auch so

% nrmse als ausschlusskriterium? aber ist vllt nicht so gut, wenn Daten
% Mist sind...

% beatwise müsste ja nur in beatwise filtered und beatwise raw vorher
% aufgeteilt werden, dann wäre das für das script hier eig kein Problem

clear all
close all
clc

%% Setup
dataset = 'CPT'; % specify dataset used
mode = 'ensemble'; % ensemble or beatwise
analyzeRef = true; % flag that decides whether or not reference is analyzed
alpha = 0.05; % significance level
activateExclusion = false; % exclude subjects based on 4 criteria
minimumCorrelation = 0.6; % minimum correlation between two ROIs needed
myFontType = 'Times New Roman'; % set font type of nice plots
myFontSize = 9; % set font size of nice plots
myFigureSize = [5 5]; % set figure size of nice plots
forceFriedman = false; % sets the statistical test used to friedman independently of number of samples
noOrthogonalContrasts = true; % flag that decides whether or not examined contrasts are orthogonal

%% Check set parameters
% mode
if(~(strcmp(mode,'ensemble')||strcmp(mode,'beatwise')))
    errordlg('Unknown mode selected. Select either ensemble or beatwise.','Input Error','modal');
    return;
end
% analyzeRef
if(~islogical(analyzeRef))
    errordlg('analyzeRef must be either true or false.','Input Error','modal');
    return;
end
% alpha
if(~((alpha==0.05) || (alpha==0.01) || (alpha==0.001)))
    errordlg('alpha needs to be at 0.05, 0.01 or 0.001.','Input Error','modal');
    return;
end
% activateExclusion
if(~islogical(activateExclusion))
    errordlg('activateExclusion must be either true or false.','Input Error','modal');
    return;
end
% minimumCorrelation
if(~(isscalar(minimumCorrelation)) && (minimumCorrelation<=1))
    errordlg('minimumCorrelation must be a scalar less than or equal to 1.','Input Error','modal');
    return;
end
% myFontType
if(~ismember(myFontType,listfonts))
    errordlg('myFontType needs to be a MATLAB font type.','Input Error','modal');
    return;
end
% myFontSize
if(~isscalar(myFontSize))
    errordlg('myFontType needs to be a MATLAB font type.','Input Error','modal');
    return;
end
% myFigureSize
if(~(size(myFigureSize,1)==1 && size(myFigureSize,2)==2 && isnumeric(myFigureSize) && isreal(myFigureSize)))
    errordlg('myFigureSize needs to be a real numeric vector of size 1x2','Input Error','modal');
    return;
end
% forceFriedman
if(~islogical(forceFriedman))
    errordlg('forceFriedman must be either true or false.','Input Error','modal');
    return;
end
% noOrthogonalContrasts
if(~islogical(noOrthogonalContrasts))
    errordlg('noOrthogonalContrasts must be either true or false.','Input Error','modal');
    return;
end

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% add path to parameter functions and cell with parameter names
addpath('..\Parameters');
% load file containing algorithms here
load('algorithmsBEST.mat','algorithms');
% check that there is an algorithm for every entry in the list
for actualAlgorithm = 1:size(algorithms,1)
    if(~(isfile(['..\Algorithms\' 'decompose' algorithms{actualAlgorithm} '.m'])))
        errordlg('Unknown decomposition algorithm selected','Input Error','modal');
        return;
    end
end
% specify folders
sourceFolder = ['..\Decomposition\Data\' dataset '\realData\'];
dataFolder = ['..\Decomposition\Data\' dataset '\decompositionParameters_' mode '\'];
resultsFolder = ['statisticalResults\' dataset '\' mode '\'];
figureFolder = ['statisticalFigures\' dataset '\' mode '\'];
% load data information
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'perfusionParametersTable.mat']);
load([sourceFolder 'epochs.mat']);
load([sourceFolder 'roiTable_RGB.mat']);
load([dataFolder 'parameterTables.mat']);
load(['..\Decomposition\Data\' dataset '\imageDummy.mat']) % load an example figure to get the size of images
patients = physiologicalMeasuresTable.SubjectID; % loads list with patients
% create new iteration of folders if there are already some existing
if(exist(figureFolder, 'dir')==7)%check if result directory exists
    existing = true;
    figureFolder = figureFolder(1:end-1);
    i = 2;
    appendix = ['_V' num2str(i) '\'];
    figureFolder = [figureFolder appendix];
    while(existing)
        if(exist(figureFolder, 'dir')~=7)
            existing = false;
        else
            figureFolder = erase(figureFolder,appendix);
            i = i+1;
            appendix = ['_V' num2str(i) '\'];
            figureFolder = [figureFolder appendix];
        end
    end
end
if(exist(resultsFolder, 'dir')==7)%check if result directory exists
    existing = true;
    resultsFolder = resultsFolder(1:end-1);
    i = 2;
    appendix = ['_V' num2str(i) '\'];
    resultsFolder = [resultsFolder appendix];
    while(existing)
        if(exist(resultsFolder, 'dir')~=7)
            existing = false;
        else
            resultsFolder = erase(resultsFolder,appendix);
            i = i+1;
            appendix = ['_V' num2str(i) '\'];
            resultsFolder = [resultsFolder appendix];
        end
    end
end

%% Set analysed epochs and variables
% set epochs which are analysed
epochsOfInterest = {'Baseline_1' 'Baseline_2' 'Baseline_3' 'AfterCPT' 'HighestSBP'}; % define here the epochs that are analysed
epochsForNorm = {'Baseline_1' 'Baseline_2' 'Baseline_3'}; % define here the epochs to which the results are normed; should be a subset of epochsOfInterest; is used for plotting
epochsForComparison = {'Baseline_3' 'AfterCPT' 'HighestSBP'}; % define here the epochs that are used for statistical assessment; should be a subset of epochsOfInterest;

% test whether epochs for comparison and normalisation are subsets of
% analyszed epochs
if(~all(ismember(epochsForNorm,epochsOfInterest)))
    errordlg('epochsForNorm are not a subset of epochsOfInterest','Setup Error','modal');
    return;
end
if(~all(ismember(epochsForComparison,epochsOfInterest)))
    errordlg('epochsForComparison are not a subset of epochsOfInterest','Setup Error','modal');
    return;
end

% set variables which are analysed
variablesOfInterest = {'b_a';'T_sys_diaMode'}; % define here the perfusion variable which should be analysed
if(analyzeRef)
    referenceParameters = {'_PP';'_rr'};
end
variablesForSortOut = {'b_a';'T_sys_diaMode'}; % here you can specify the variables which are considered to sort out epochs (should contain variablesOfInterest but can contain more)

% test whether variables for sort out contain analysed variables
if(~all(ismember(variablesOfInterest,variablesForSortOut)))
    errordlg('variablesForSortOut does not contain all variablesOfInterest','Setup Error','modal');
    return;
end

%% Initialize variables to be saved
combinedTables = struct; % initialize struct for storing combined tables
physMeadTables = struct; % initialize struct for storing physiological measurement tables

% initialize matrix of p values
pValue = NaN(size(variablesOfInterest,1),size(algorithms,1));
pValueCorrected = NaN(size(variablesOfInterest,1),size(algorithms,1));
hGes = NaN(size(variablesOfInterest,1),size(algorithms,1));
if(analyzeRef)
    % initalize matrix of p values (reference parameters)
    pValueRef = NaN(size(referenceParameters,1),size(algorithms,1));
    pValueTotalCorrected = NaN(size(variablesOfInterest,1)+size(referenceParameters,1),size(algorithms,1));
    hGes = NaN(size(variablesOfInterest,1)+size(referenceParameters,1),size(algorithms,1));
end

% initialize exclusion criteria
exclusionCriteria = zeros(size(algorithms,1),4);
removeIndices = cell(size(algorithms,1),1);

% create cell array to store the tables of p values of the post-hoc tests
% of all variables of interest and algorithms
pTables = cell(size(variablesOfInterest,1),size(algorithms,1));
if(analyzeRef)
    % create cell array to store the tables of p values of the post-hoc tests
    % of all variables of interest and algorithms
    pTablesRef = cell(size(referenceParameters,1),size(algorithms,1));
end

for actualAlgorithm = 1:size(algorithms,1)
    %% Create possbily missing paths
    % figure folders
    if(exist([figureFolder algorithms{actualAlgorithm} '\'], 'dir')~=7)%check if result directory exists
        mkdir([figureFolder algorithms{actualAlgorithm} '\'])%if not, create it
    end
    
    % results folder
    if(exist(resultsFolder, 'dir')~=7)%check if result directory exists
        mkdir(resultsFolder)%if not, create it
    end
    
    %% Shorten tables to relevant information
    % combine perfusion parameter table and decomposition parameter table
    combinedTable = join(perfusionParametersTable,parameterTables.(algorithms{actualAlgorithm}));
    columnNames = combinedTable.Properties.VariableNames;
    desiredOrder = columnNames(1);
    for currentEpoch = 1:numel(epochs)
        nextEntries = columnNames(contains(columnNames,epochs{currentEpoch,1}));
        desiredOrder = [desiredOrder nextEntries];
    end
    combinedTable = combinedTable(:,desiredOrder);
    
    % make temporal copies of tables
    physiologicalMeasuresTableTmp = physiologicalMeasuresTable; % temporal copy of reference parameters
    combinedTableTmp = combinedTable; % temporal copy of perfusion parameters
    
    % find columns in the perfusion parameters for analysis (predefined epoch and variable)
    % each variable in different row in order to be able to distinguish
    % variables
    columnsOfInterestPerfusion = false(size(variablesOfInterest,1),size(combinedTableTmp,2)); % preallocate memory for columnsOfInterestPerfusion
    for currentVariableOfInterest = 1:size(variablesOfInterest,1)
        columnsOfInterestPerfusion(currentVariableOfInterest,:) = ...
            endsWith(combinedTableTmp.Properties.VariableNames,variablesOfInterest(currentVariableOfInterest)) & ...
            contains(combinedTableTmp.Properties.VariableNames,epochsOfInterest); % also muss Variablennamen + Epoche beinhalten --> meine Tables haben gleiches Format
    end
    
    % find columns in the reference parameters for analysis (predefined epoch and variable)
    % each variable in different row in order to be able to distinguish
    % variables
    if(analyzeRef)
        columnsOfInterestReference = false(size(referenceParameters,1),size(physiologicalMeasuresTableTmp,2)); % preallocate memory for columnsOfInterestReference
        for currentReferenceParameter = 1:size(referenceParameters,1)
            columnsOfInterestReference(currentReferenceParameter,:) =...
                endsWith(physiologicalMeasuresTableTmp.Properties.VariableNames,referenceParameters(currentReferenceParameter)) & ...
                contains(physiologicalMeasuresTableTmp.Properties.VariableNames,epochsOfInterest);
        end
    end
    
    % define columns in the perfusion parameters for sort out (predefined epoch and variable)
    columnsOfInterestSortOut = ...
        endsWith(combinedTableTmp.Properties.VariableNames,variablesForSortOut) & ...
        contains(combinedTableTmp.Properties.VariableNames,epochsOfInterest);
    % define columns in the physiologiccal measures for sort out (predefined epoch and variable)
    columnsOfInterestRefSortOut = ...
        endsWith(physiologicalMeasuresTableTmp.Properties.VariableNames,variablesForSortOut) & ...
        contains(physiologicalMeasuresTableTmp.Properties.VariableNames,epochsOfInterest);
    
    %% Exclusion criteria
    % list exclusion criteria here:
    % (1):  placement of analysis intervals failed or no stable forehead region
    %       could be defined for strong subject movements
    % (2):  template generation failed due to too few beats in the iPPG signal
    %       or in the reference PPG (for beatwise analysis: too few beats for
    %       calculation of mean)
    % (3):  strong displacements of the subjects’ heads between analysis
    %       intervals
    % (4):  ROI contained less than 1000 pixels
    
    % criterium 1
    roiTableMissing = ismissing(roiTable,{NaN});
    rowsWithMissing = roiTable(any(roiTableMissing,2),:);
    exclusionCriteria(actualAlgorithm,1) = size(rowsWithMissing,1); % subjects that are excluded due to (1)
    
    % identify rows that change ROIs between epochs significantly
    excludeSubjects = [];
    roiMaskSize = zeros(size(roiTable,1),numel(epochsOfInterest)); % initialize matrix
    roiMaskSize(roiMaskSize == 0) = NaN; % initialize matrix with NaNs
    for actualPatientNumber = 1:size(roiTable,1) % go thorugh all subjects
        exitFlag = 0; % set back exitFlag for new subject
        roiMask = []; % set back roiMask for new subject
        correlation = []; % set back roiMask for new subject
        for currentInterval = 1:numel(epochsOfInterest) % get ROI mask for the intervals
            roiPoly = roiTable.([epochsOfInterest{currentInterval} '_ROI']){actualPatientNumber}; % get polygon of ROI
            if(isnan(roiPoly)) % skip the subject if no ROI is available
                exitFlag = 1;
                excludeSubjects(end+1) = actualPatientNumber; % exclude subjects that have a missing ROI
                break
            end
            roiMask(:,:,currentInterval) = poly2mask(roiPoly(:,1),roiPoly(:,2),size(imageDummy,1),size(imageDummy,2)); % get ROI mask
        end
        for currentRoiMask = 1:size(roiMask,3)
            roiMaskSize(actualPatientNumber,currentRoiMask) = sum(sum(roiMask(:,:,currentRoiMask))); % get size of ROI
        end
        if(exitFlag==1)
            continue
        end
        for i = 1:currentInterval % calculate correlation between roi maskes
            for j = 1:currentInterval
                correlation(i,j) = corr2(roiMask(:,:,i),roiMask(:,:,j));
            end
        end
        if((min(min(correlation))<minimumCorrelation) | (~isempty(find(isnan(correlation)))))
            excludeSubjects(end+1)=actualPatientNumber;
        end
    end
    
    % find subjects to be removed and number of subjects for each exclusion criterium
    [NaNrows,NaNcols] = find(isnan(combinedTableTmp{:,columnsOfInterestSortOut})); % find rows with NaNs - here only the variable of interest or more variables can be used
    NaNrows = unique(NaNrows);
    if(analyzeRef)
        [NaNrowsRef,NaNcolsRef] = find(isnan(physiologicalMeasuresTableTmp{:,columnsOfInterestRefSortOut}));
        NaNrowsRef = unique(NaNrowsRef);
        removeIndices{actualAlgorithm,1} = unique([excludeSubjects (find(any(roiMaskSize<1000,2)))' NaNrows' NaNrowsRef']); % get indices to remove
        exclusionCriteria(actualAlgorithm,2) = numel(unique([NaNrows; NaNrowsRef])) - size(rowsWithMissing,1); % subjects that are excluded due to (2)
        exclusionCriteria(actualAlgorithm,3) = numel(excludeSubjects); % subjects that are excluded due to (3)
        exclusionCriteria(actualAlgorithm,4) = numel(find(any(roiMaskSize<1000,2))); % subjects that are excluded due to (4)
    else
        removeIndices{actualAlgorithm,1} = unique([excludeSubjects (find(any(roiMaskSize<1000,2)))' NaNrows']); % get indices to remove
        exclusionCriteria(actualAlgorithm,2) = numel(NaNrows) - size(rowsWithMissing,1); % subjects that are excluded due to (2)
        exclusionCriteria(actualAlgorithm,3) = numel(excludeSubjects); % subjects that are excluded due to (3)
        exclusionCriteria(actualAlgorithm,4) = numel(find(any(roiMaskSize<1000,2))); % subjects that are excluded due to (4)
    end
    
    % remove subjects from tables & store final tables in struct
    if(activateExclusion == true)
        combinedTableTmp(removeIndices{actualAlgorithm,1},:) = []; % remove subjects that are sorted out due to exclusion criteria
        physiologicalMeasuresTableTmp(removeIndices{actualAlgorithm,1},:) = []; % remove subjects that are sorted out due to exclusion criteria
    else
        if(analyzeRef)
            removeIndices{actualAlgorithm,1} = unique([NaNrows' NaNrowsRef']); % only exclude subjects with missing values
        else
            removeIndices{actualAlgorithm,1} = NaNrows'; % only exclude subjects with missing values
        end
        combinedTableTmp(removeIndices{actualAlgorithm,1},:) = []; % remove subjects that are sorted out due to exclusion criteria
        physiologicalMeasuresTableTmp(removeIndices{actualAlgorithm,1},:) = []; % remove subjects that are sorted out due to exclusion criteria
    end
    combinedTables.(algorithms{actualAlgorithm}) = combinedTableTmp;
    physMeasTables.(algorithms{actualAlgorithm}) = physiologicalMeasuresTableTmp;
    if(size(combinedTableTmp,1)<2)
        continue % skip statistical calculation if friedman test is not possible for this algorithm (2 rows are mandatory)
    end
    
    %% do omnibus analysis
    % initialize variableNames --> kann rausgekürzt werden?
    if(analyzeRef)
        variableNames = cell(size(variablesOfInterest,1)+size(referenceParameters,1),1);
    else
        variableNames = cell(size(variablesOfInterest,1),1);
    end
    
    %% Do omnibus analysis of camera parameters
    % define column which are used to normalize the result patient wise
    columnsOfInterestNorm = false(size(variablesOfInterest,1),size(combinedTableTmp,2));
    normValue = zeros(size(combinedTableTmp,1),size(variablesOfInterest,1));
    parameterNorm = zeros(size(combinedTableTmp,1),size(epochsOfInterest,2),size(variablesOfInterest,1));
    for currentVariableOfInterest = 1:size(variablesOfInterest,1)
        columnsOfInterestNorm(currentVariableOfInterest,:)=...
            endsWith(combinedTableTmp.Properties.VariableNames,variablesOfInterest(currentVariableOfInterest)) & ...
            contains(combinedTableTmp.Properties.VariableNames,epochsForNorm);
        normValue(:,currentVariableOfInterest) = ...
            mean(combinedTableTmp{:,columnsOfInterestNorm(currentVariableOfInterest,:)},2); % calculate value to which variable is normalized
        parameterNorm(:,:,currentVariableOfInterest) = ...
            combinedTableTmp{:,columnsOfInterestPerfusion(currentVariableOfInterest,:)}./repmat(normValue(:,currentVariableOfInterest),[1 size(combinedTableTmp{:,columnsOfInterestPerfusion(currentVariableOfInterest,:)},2)]); % normalize variable
        currentVariableName = erase(variablesOfInterest{currentVariableOfInterest},'_'); % delete '_' from variable name to avoid first char to become subscript
        variableNames{currentVariableOfInterest,1} = join([currentVariableName,'_iPPG']); % save all better Strings as a whole
        
        %boxplot
        figure;
        boxplot(parameterNorm(:,:,currentVariableOfInterest))
        ylim auto
        xticks(1:size(epochsOfInterest,2))
        set(gca,'TickLabelInterpreter', 'tex'); % '_' is interpreted as subscript
        xticklabels(epochsOfInterest)
        xtickangle(45)
        ylabel([currentVariableName '/a.u.'])
        title([currentVariableName ' ' algorithms{actualAlgorithm} ' iPPG']);
        savefig([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_iPPG_Box.fig'])
        close
        
        % nice boxplot
        figure;
        boxplot(parameterNorm(:,:,currentVariableOfInterest),'Colors','k','Widths',0.3)
        h = findobj(gca,'tag','Outliers');
        delete(h) % delete outliers from plot
        ylim auto
        set(gca,'TickLabelInterpreter', 'tex'); % '_' is interpreted as subscript
        xticklabels(epochsOfInterest)
        xtickangle(45)
        ylabel([currentVariableName '/a.u.'])
        title([currentVariableName ' ' algorithms{actualAlgorithm} ' iPPG']);
        box off
        set(gcf, 'Units', 'centimeters');
        set(gcf, 'PaperUnits', 'centimeters');
        currentPos=get(gcf,'Position');
        set(gcf, 'Position', [currentPos(1) currentPos(2) myFigureSize]);
        set(findobj(gcf,'type','axes'),...
            'FontSize', myFontSize,...
            'FontName', myFontType,...
            'FontWeight','normal',...
            'TitleFontWeight','normal');
        savefig([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_iPPG_niceBox.fig'])
        close
        
        % plot
        figure;
        plot(parameterNorm(:,:,currentVariableOfInterest)')
        ylim auto
        xticks(1:size(epochsOfInterest,2))
        xticklabels(epochsOfInterest)
        xtickangle(45)
        ylabel([currentVariableName '/a.u.'])
        title([currentVariableName ' ' algorithms{actualAlgorithm} ' iPPG']);
        savefig([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_iPPG_Line.fig'])
        close
        
        if((size(parameterNorm,1)>=25) && ~(forceFriedman))
            % repeated measures anova
            anovaTable = array2table(parameterNorm(:,:,currentVariableOfInterest),'VariableNames',epochsOfInterest); % create table from input matrix
            wilkNot = []; % initialize Wilkinson notation
            for entry = 1:size(epochsOfInterest,2)
                wilkNot = [wilkNot epochsOfInterest{1,entry} ',']; % specify Wilkinson notation
            end
            wilkNot(end) = ''; % delete last comma
            wilkNot = [wilkNot ' ~ 1']; % specify that there is no predictor variable in repeated measures model
            %withinDesign = cell2table(epochsOfInterest'); % specify within design subject factors
            %rm = fitrm(anovaTable,wilkNot,'WithinDesign',withinDesign); % fit repeated measures model
            rm = fitrm(anovaTable,wilkNot); % fit repeated measures model
            ranovatbl = ranova(rm); % perform repeated measures anova
            pValue(currentVariableOfInterest,actualAlgorithm) = ranovatbl.pValue(1);
        else
            % friedman test
            pValue(currentVariableOfInterest,actualAlgorithm) = friedman(parameterNorm(:,:,currentVariableOfInterest),1,'off');
        end
    end
    %% If required, do omnibus analysis of reference parameters
    if(analyzeRef)
        % define column which are used to normalize the result patient wise
        columnsOfInterestNormReference = false(size(referenceParameters,1),size(physiologicalMeasuresTableTmp,2));
        normValueRef = zeros(size(physiologicalMeasuresTableTmp,1),size(referenceParameters,1));
        parameterNormRef = zeros(size(physiologicalMeasuresTableTmp,1),size(epochsOfInterest,2),size(referenceParameters,1));
        for currentReferenceParameter = 1:size(referenceParameters,1)
            columnsOfInterestNormReference(currentReferenceParameter,:)=...
                endsWith(physiologicalMeasuresTableTmp.Properties.VariableNames,referenceParameters(currentReferenceParameter)) & ...
                contains(physiologicalMeasuresTableTmp.Properties.VariableNames,epochsForNorm);
            normValueRef(:,currentReferenceParameter)=mean(physiologicalMeasuresTableTmp{:,columnsOfInterestNormReference(currentReferenceParameter,:)},2);
            parameterNormRef(:,:,currentReferenceParameter)=physiologicalMeasuresTableTmp{:,columnsOfInterestReference(currentReferenceParameter,:)}./repmat(normValueRef(:,currentReferenceParameter),[1 size(physiologicalMeasuresTableTmp{:,columnsOfInterestReference(currentReferenceParameter,:)},2)]);
            currentVariableName = erase(referenceParameters{currentReferenceParameter},'_'); % delete '_' from variblae name to avoid first char to become subscript
            variableNames{size(variablesOfInterest,1)+currentReferenceParameter,1} = join([currentVariableName,'_ref']); % save all better Strings as a whole
            
            % boxplot
            figure;
            boxplot(parameterNormRef(:,:,currentReferenceParameter))
            ylim auto
            xticks(1:size(epochsOfInterest,2))
            set(gca,'TickLabelInterpreter', 'tex'); % '_' is interpreted as subscript
            xticklabels(epochsOfInterest)
            xtickangle(45)
            ylabel([currentVariableName '/a.u.'])
            title([currentVariableName ' ' algorithms{actualAlgorithm}])
            savefig([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_ref_Box.fig'])
            close
            
            % nice boxplot
            figure;
            boxplot(parameterNormRef(:,:,currentReferenceParameter),'Colors','k','Widths',0.3)
            h = findobj(gca,'tag','Outliers');
            delete(h) % delete outliers from plot
            ylim auto
            set(gca,'TickLabelInterpreter', 'tex'); % '_' is interpreted as subscript
            xticklabels(epochsOfInterest)
            xtickangle(45)
            ylabel([currentVariableName '/a.u.'])
            title([currentVariableName ' ' algorithms{actualAlgorithm}]);
            box off
            set(gcf, 'Units', 'centimeters');
            set(gcf, 'PaperUnits', 'centimeters');
            currentPos=get(gcf,'Position');
            set(gcf, 'Position', [currentPos(1) currentPos(2) myFigureSize]);
            set(findobj(gcf,'type','axes'),...
                'FontSize', myFontSize,...
                'FontName', myFontType,...
                'FontWeight','normal',...
                'TitleFontWeight','normal');
            savefig([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_ref_niceBox.fig'])
            close
            
            % plot
            figure;
            plot(parameterNormRef(:,:,currentReferenceParameter)')
            ylim auto
            xticks(1:size(epochsOfInterest,2))
            xticklabels(epochsOfInterest)
            xtickangle(45)
            ylabel([currentVariableName '/a.u.'])
            title([currentVariableName ' ' algorithms{actualAlgorithm}])
            savefig([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_ref_Line.fig'])
            close
            
            
            if((size(parameterNormRef,1)>=25) && ~(forceFriedman))
                % repeated measures anova
                anovaTable = array2table(parameterNormRef(:,:,currentReferenceParameter),'VariableNames',epochsOfInterest); % create table from input matrix
                wilkNot = []; % initialize Wilkinson notation
                for entry = 1:size(epochsOfInterest,2)
                    wilkNot = [wilkNot epochsOfInterest{1,entry} ',']; % specify Wilkinson notation
                end
                wilkNot(end) = ''; % delete last comma
                wilkNot = [wilkNot ' ~ 1']; % specify that there is no predictor variable in repeated measures model
                %withinDesign = cell2table(epochsOfInterest'); % specify within design subject factors
                %rm = fitrm(anovaTable,wilkNot,'WithinDesign',withinDesign); % fit repeated measures model
                rm = fitrm(anovaTable,wilkNot); % fit repeated measures model
                ranovatbl = ranova(rm); % perform repeated measures anova
                pValueRef(currentReferenceParameter,actualAlgorithm) = ranovatbl.pValue(1);
            else
                % friedman test
                pValueRef(currentReferenceParameter,actualAlgorithm) = friedman(parameterNormRef(:,:,currentReferenceParameter),1,'off');
            end
        end
        
        % display original and corrected p values of Friedman test
        pValueTotal = [pValue(:,actualAlgorithm); pValueRef(:,actualAlgorithm)];
        [pValueTotalCorrected(:,actualAlgorithm), hGes(:,actualAlgorithm)] = bonf_holm(pValueTotal,alpha);
        disp('##################');
        disp([algorithms{actualAlgorithm} ':']);
        disp('##################');
        for currentVariableOfInterest = 1:size(variablesOfInterest,1)
            disp(variableNames{currentVariableOfInterest,1});
            disp(['p = ' num2str(pValue(currentVariableOfInterest,actualAlgorithm))]);
            disp(['pCorrected = ' num2str(pValueTotalCorrected(currentVariableOfInterest,actualAlgorithm))]);
            fprintf('\n'); % print empty line for overview purposes
        end
        for currentReferenceParameter = 1:size(referenceParameters,1)
            disp(variableNames{size(variablesOfInterest,1)+currentReferenceParameter,1});
            disp(['p = ' num2str(pValueRef(currentReferenceParameter,actualAlgorithm))]);
            disp(['pCorrected = ' num2str(pValueTotalCorrected(size(variablesOfInterest,1)+currentReferenceParameter,actualAlgorithm))]);
            fprintf('\n'); % print empty line for overview purposes
        end
    else
        [pValueCorrected(:,actualAlgorithm), hGes(:,actualAlgorithm)] = bonf_holm(pValue(:,actualAlgorithm),alpha);
        disp('##################');
        disp([algorithms{actualAlgorithm} ':']);
        disp('##################');
        for currentVariableOfInterest = 1:size(variablesOfInterest,1)
            disp(variableNames{currentVariableOfInterest,1});
            disp(['p = ' num2str(pValue(currentVariableOfInterest,actualAlgorithm))]);
            disp(['pCorrected = ' num2str(pValueCorrected(currentVariableOfInterest,actualAlgorithm))]);
            fprintf('\n'); % print empty line for overview purposes
        end
    end
    
    %% Do post-hoc analysis of perfusion parameters
    for currentVariableOfInterest = 1:size(variablesOfInterest,1)
        pMatrix = zeros(numel(epochsForComparison));
        if(hGes(currentVariableOfInterest,actualAlgorithm))  % use h to determine if p < 0.05
            for i=1:numel(epochsForComparison)
                for j=1:numel(epochsForComparison)
                    if((size(parameterNorm,1)>=25) && ~(forceFriedman))
                        [~,pMatrix(i,j)] = ttest2(parameterNorm(:, ...
                            find(contains(epochsOfInterest,epochsForComparison{i})), ...
                            currentVariableOfInterest),...
                            parameterNorm(:,find(contains(epochsOfInterest,epochsForComparison{j})), ...
                            currentVariableOfInterest));
                    else
                        pMatrix(i,j) = signrank(parameterNorm(:, ...
                            find(contains(epochsOfInterest,epochsForComparison{i})), ...
                            currentVariableOfInterest),...
                            parameterNorm(:,find(contains(epochsOfInterest,epochsForComparison{j})), ...
                            currentVariableOfInterest));
                    end
                end
            end
            if(noOrthogonalContrasts)
                pPostHoc = zeros(nchoosek(numel(epochsForComparison),2),1);
                counter = 1;
                for i=1:numel(epochsForComparison)
                    for j=1:numel(epochsForComparison)
                        if(i < j)
                            pPostHoc(counter) = pMatrix(i,j);
                            counter = counter + 1;
                        end
                    end
                end
                [pPostHocCorrected,~] = bonf_holm(pPostHoc,alpha);
                for currentPValue = 1:numel(pPostHoc)
                    pMatrix(pMatrix==pPostHoc(currentPValue)) = pPostHocCorrected(currentPValue);
                end
            end
        else
            for i=1:numel(epochsForComparison)
                for j=1:numel(epochsForComparison)
                    pMatrix(i,j) = NaN;
                end
            end
        end
        pTable = table('Size',[numel(epochsForComparison) numel(epochsForComparison)],...
            'VariableTypes',repmat({'double'},[1 numel(epochsForComparison)]),...
            'VariableNames',epochsForComparison,...
            'RowNames',epochsForComparison);
        pTable.Properties.Description = ...
            ['p values for ' variablesOfInterest{currentVariableOfInterest} ' of iPPG'];
        for i=1:numel(epochsForComparison)
            for j=1:numel(epochsForComparison)
                pTable{i,j} = pMatrix(i,j);
            end
        end
        pTables{currentVariableOfInterest,actualAlgorithm} = pTable;
    end
    
    %% If required, do post-hoc analysis of reference parameters
    if(analyzeRef)
        for currentReferenceParameter = 1:size(referenceParameters,1)
            pMatrixRef = zeros(numel(epochsForComparison));
            if(hGes(size(variablesOfInterest,1)+currentReferenceParameter,actualAlgorithm))  % use h to determine if p < 0.05
                for i=1:numel(epochsForComparison)
                    for j=1:numel(epochsForComparison)
                        if((size(parameterNormRef,1)>=25) && ~(forceFriedman))
                            [~,pMatrixRef(i,j)] = ttest2(parameterNormRef(:, ...
                            find(contains(epochsOfInterest,epochsForComparison{i})), ...
                            currentReferenceParameter), ...
                            parameterNormRef(:,find(contains(epochsOfInterest,epochsForComparison{j})), ...
                            currentReferenceParameter));
                        else
                            pMatrixRef(i,j) = signrank(parameterNormRef(:, ...
                            find(contains(epochsOfInterest,epochsForComparison{i})), ...
                            currentReferenceParameter), ...
                            parameterNormRef(:,find(contains(epochsOfInterest,epochsForComparison{j})), ...
                            currentReferenceParameter));
                        end
                    end
                end
                if(noOrthogonalContrasts)
                    pPostHoc = zeros(nchoosek(numel(epochsForComparison),2),1);
                    counter = 1;
                    for i=1:numel(epochsForComparison)
                        for j=1:numel(epochsForComparison)
                            if(i < j)
                                pPostHoc(counter) = pMatrixRef(i,j);
                                counter = counter + 1;
                            end
                        end
                    end
                    [pPostHocCorrected,~] = bonf_holm(pPostHoc,alpha);
                    for currentPValue = 1:numel(pPostHoc)
                        pMatrixRef(pMatrixRef==pPostHoc(currentPValue)) = pPostHocCorrected(currentPValue);
                    end
                end
            else
                for i=1:numel(epochsForComparison)
                    for j=1:numel(epochsForComparison)
                        pMatrixRef(i,j) = NaN;
                    end
                end
            end
            pTableRef = table('Size',[numel(epochsForComparison) numel(epochsForComparison)],...
                'VariableTypes',repmat({'double'},[1 numel(epochsForComparison)]),...
                'VariableNames',epochsForComparison,...
                'RowNames',epochsForComparison);
            pTableRef.Properties.Description = ...
                ['p values for ' referenceParameters{currentReferenceParameter} ' of Reference PPG'];
            for i=1:numel(epochsForComparison)
                for j=1:numel(epochsForComparison)
                    pTableRef{i,j} = pMatrixRef(i,j);
                end
            end
            pTablesRef{currentReferenceParameter,actualAlgorithm} = pTableRef;
        end
    end
end

%% Save statistical results
% create table store number of excluded subjects as a whole and for every criterium
exclusionCriteria = array2table(exclusionCriteria,...
    'RowNames',algorithms,...
    'VariableNames',{'Criterion1','Criterion2','Criterion3','Criterion4'});
save([resultsFolder 'exclusionCriteria.mat'],...
    'exclusionCriteria','removeIndices')

% save tables with parameters after exclusion
save([resultsFolder 'adjustedParamTables.mat'],'combinedTables','physMeasTables')

% create table and store omnibus test results
if(analyzeRef)
    pValue = array2table(pValue,'RowNames',variablesOfInterest,'VariableNames',algorithms);
    pValueRef = array2table(pValueRef,'RowNames',referenceParameters,'VariableNames',algorithms);
    pValueTotalCorrected = array2table(pValueTotalCorrected,'RowNames',variableNames,'VariableNames',algorithms);
    save([resultsFolder '\omnibusTest.mat'],...
        'pValue','pValueTotalCorrected','pValueRef','hGes')
else
    pValue = array2table(pValue,'RowNames',variableNames,'VariableNames',algorithms);
    pValueCorrected = array2table(pValueCorrected,'RowNames',variableNames,'VariableNames',algorithms);
    save([resultsFolder '\omnibusTest.mat'],...
        'pValue','pValueCorrected','hGes')
end

% create table and save post-hoc p values
pTables = cell2table(pTables,'VariableNames',algorithms);
save([resultsFolder 'postHociPPG.mat'],'pTables')
if(analyzeRef)
    pTablesRef = cell2table(pTablesRef,'VariableNames',algorithms);
    save([resultsFolder 'postHocRef.mat'],'pTablesRef')
end

%% Put significance bars into box plots
% iPPG parameters
for actualAlgorithm = 1:size(algorithms,1)
    for currentVariableOfInterest = 1:size(variablesOfInterest,1)
        currentVariableName = variableNames{currentVariableOfInterest,1};
        if(hGes(currentVariableOfInterest,actualAlgorithm)==1)
            [row,col] = find(pTables.(algorithms{actualAlgorithm,1}){currentVariableOfInterest,1}{:,:}<alpha);
            if(~isempty(row))
                indices = unique(sort([row col],2),'rows');
                values = zeros(1,size(indices,1));
                for currentIndices = 1:size(indices,1)
                    values(1,currentIndices) = pTables.(algorithms{actualAlgorithm,1}){currentVariableOfInterest,1}{indices(currentIndices,1),indices(currentIndices,2)};
                    indices(currentIndices,1) = find(strcmp(epochsOfInterest,epochsForComparison{indices(currentIndices,1)}));
                    indices(currentIndices,2) = find(strcmp(epochsOfInterest,epochsForComparison{indices(currentIndices,2)}));
                end
                pairs = (num2cell(indices,2))';
                currentFigure = [figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_Box.fig'];
                openfig(currentFigure);
                hold on
                sigstar(pairs,values);
                hold off
                savefig(currentFigure);
                matlab2tikz([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_Box.tex'])
                close
                currentFigure = [figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_niceBox.fig'];
                openfig(currentFigure);
                hold on
                sigstar(pairs,values);
                hold off
                savefig(currentFigure);
                matlab2tikz([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_niceBox.tex'])
                close
            end
        end
    end
end

% reference parameters
if(analyzeRef)
    for actualAlgorithm = 1:size(algorithms,1)
        for currentReferenceParameter = 1:size(referenceParameters,1)
            currentVariableName = variableNames{size(variablesOfInterest,1)+currentReferenceParameter,1};
            if(hGes(size(variablesOfInterest,1)+currentReferenceParameter,actualAlgorithm)==1)
                [row,col] = find(pTablesRef.(algorithms{actualAlgorithm,1}){currentReferenceParameter,1}{:,:}<alpha);
                if(~isempty(row))
                    indices = unique(sort([row col],2),'rows');
                    values = zeros(1,size(indices,1));
                    for currentIndices = 1:size(indices,1)
                        values(1,currentIndices) = pTablesRef.(algorithms{actualAlgorithm,1}){currentReferenceParameter,1}{indices(currentIndices,1),indices(currentIndices,2)};
                        indices(currentIndices,1) = find(strcmp(epochsOfInterest,epochsForComparison{indices(currentIndices,1)}));
                        indices(currentIndices,2) = find(strcmp(epochsOfInterest,epochsForComparison{indices(currentIndices,2)}));
                    end
                    pairs = (num2cell(indices,2))';
                    currentFigure = [figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_Box.fig'];
                    openfig(currentFigure);
                    hold on
                    sigstar(pairs,values);
                    hold off
                    savefig(currentFigure);
                    matlab2tikz([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_Box.tex'])
                    close
                    currentFigure = [figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_niceBox.fig'];
                    openfig(currentFigure);
                    hold on
                    sigstar(pairs,values);
                    hold off
                    savefig(currentFigure);
                    matlab2tikz([figureFolder algorithms{actualAlgorithm} '\' currentVariableName '_niceBox.tex'])
                    close
                end
            end
        end
    end
end