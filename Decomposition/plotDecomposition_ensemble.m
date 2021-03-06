clear all
close all
clc

dataset='CPT';
sourceFolder=['Data\' dataset '\realData\'];
dataFolder=['Data\' dataset '\results_ensemble\'];
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']);
load('algorithmsBESTTEST.mat');
patients=physiologicalMeasuresTable.SubjectID;


for actualPatientNumber=1:size(patients,1)
    fileID = patients{actualPatientNumber};
    %% loop over epochs & algorithms
    for actualAlgorithm = 1:size(algorithms,1)
        %% init figure and subplots for output
        f1=figure('PaperUnits', 'normalized',...
            'PaperPosition',[0 0 1 1],...
            'NumberTitle', 'off',...
            'Name',[fileID]);
        algorithmName = algorithms{actualAlgorithm};
        for currentInterval=1:6
            currentFile = [dataFolder fileID '\' epochs{currentInterval} '\' algorithmName];
            try%try loading
                load(currentFile)
            catch%if not loadable
                meanSeg=NaN;
                y=cell(3,1);
                signal_mod=NaN;
            end
            
            %% show signal
            myaxSignal(currentInterval)=subplot(4,3,currentInterval);
            if(~isnan(meanSeg))
                plot(meanSeg);
            else
                plot(1:100,0);%create dummy plot
            end
            title({[fileID ' ' epochs{currentInterval}]; algorithmName; ['beats: ', num2str(numel(beatIndicesEnsembleBeat))]},'Interpreter','none')
            
            
            %% show decomposition
            myaxSignal(currentInterval+6)=subplot(4,3,currentInterval+6);%myaxSignal(currentInterval+6)
            if(~isempty(y{1}))
                plot(beatStartIndex:beatStopIndex,meanSeg(beatStartIndex:beatStopIndex),'LineWidth',3); hold on
                plot(beatStartIndex:beatStopIndex,signal_mod,'Color','r')
                for i = 1:size(y,1)
                    plot(beatStartIndex:beatStopIndex,y{i},'k')
                end
                hold off
            else
                plot(1:numel(meanSeg),0)
            end
            title({[fileID ' ' epochs{currentInterval}]; algorithmName; ['beats: ', num2str(numel(beatIndicesEnsembleBeat))]},'Interpreter','none')
            
            
            
        end
        %% format subplots
        set(myaxSignal,'XLim', [0 max(max(cell2mat({myaxSignal.XLim}')))])
        set(myaxSignal(1:12),'YLim',[min(min(cell2mat({myaxSignal(1:12).YLim}'))) max(max(cell2mat({myaxSignal(1:12).YLim}'))) ])
        
        %% store plot
        print([dataFolder 'decompositionPlot_' algorithmName],'-dpsc','-append')
        
        %close plot
        close all
        clear myaxSignal myaxModel f1
    end 
end
