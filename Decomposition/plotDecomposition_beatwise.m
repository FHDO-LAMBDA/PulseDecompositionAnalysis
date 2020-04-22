clear all
close all
clc
%
% shows only 6 of 7 epochs due to optimal using of a pdf page
%
% can't this be combined with the ensemble script?
%

dataset='CPT';
sourceFolder=['Data\' dataset '\realData\'];
dataFolder=['Data\' dataset '\results_beatwise\'];
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']);
load('algorithmsSORELLI.mat');
patients=physiologicalMeasuresTable.SubjectID;


for actualPatientNumber=1:size(patients,1)
    fileID = patients{actualPatientNumber};
    %% loop over epochs & beats & algorithms
    for currentInterval=1:6
        for actualAlgorithm = 1:size(algorithms,1)
            %% init figure and subplots for output
            f1=figure('PaperUnits', 'normalized',...
                'PaperPosition',[0 0 1 1],...
                'NumberTitle', 'off',...
                'Name',[fileID '_' epochs{currentInterval} '_' algorithms{actualAlgorithm,1}]);
            algorithmName = algorithms{actualAlgorithm};
            currentFile = [dataFolder fileID '\' epochs{currentInterval} '\' algorithmName];
            try%try loading
                load(currentFile)
            catch%if not loadable
                decompositionResults.singleBeats = NaN;
                decompositionResults.meanSeg = NaN;
                decompositionResults.y=cell(3,1);
                decompositionResults.signal_mod=NaN;
                beatIndicesEnsembleBeat = NaN;
            end
            % only space for 12 decompositions --> display only first 12
            % beats
            if(length(beatIndicesEnsembleBeat)>12)
                beatIndicesEnsembleBeat = beatIndicesEnsembleBeat(1:12);
            end
            
            %% show signal
            for currentBeat = 1:length(beatIndicesEnsembleBeat)
                %% show decomposition
                myaxSignal(currentBeat)=subplot(4,3,currentBeat);
                if(~isnan(decompositionResults(currentBeat).singleBeats))
                    plot(decompositionResults(currentBeat).singleBeats,'LineWidth',3); hold on
                    plot(decompositionResults(currentBeat).signal_mod,'Color','r')
                    for i = 1:size(decompositionResults(currentBeat).y,1)
                        plot(decompositionResults(currentBeat).y{i},'k')
                    end
                    hold off
                else
                    plot(1:numel(decompositionResults(currentBeat).singleBeats),0) % show dummy plot
                end
                title({[fileID ' ' epochs{currentInterval}]; algorithmName; ['beat: ', num2str(currentBeat), '/', num2str(numel(beatIndicesEnsembleBeat))]},'Interpreter','none')
            end
            clear decompositionResults
            %% format subplots
            set(myaxSignal,'XLim', [0 max(max(cell2mat({myaxSignal.XLim}')))])
            set(myaxSignal(1:length(beatIndicesEnsembleBeat)),'YLim',[min(min(cell2mat({myaxSignal(1:length(beatIndicesEnsembleBeat)).YLim}'))) max(max(cell2mat({myaxSignal(1:length(beatIndicesEnsembleBeat)).YLim}'))) ])
            
            %% store plot
            print([dataFolder 'decompositionPlot_' algorithmName],'-dpsc','-append')
            
            %close plot
            close all
            clear myaxSignal myaxModel f1
        end
    end
end
