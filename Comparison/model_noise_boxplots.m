clc
clear all
close all

% NaN is omitted

Name = 'SET5Kernel';

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% load file containing Algorithmnames here
load(['algorithms' Name '.mat'],'algorithms');
% load file specifying colors for the plots here
load('colors2alg.mat','colors');
% check that there is an algorithm for every entry in the list
for actualAlgorithm = 1:size(algorithms,1)
    if(~(isfile(['..\Algorithms\' 'decompose' algorithms{actualAlgorithm} '.m'])))
        errordlg('Unknown decomposition algorithm selected','Input Error','modal');
        return
    end
end
% check if there is a color for every algorithm
if(~(size(algorithms,1) == size(colors,1)))
    errordlg('Number of colors and algorithms specified does not match','Input Error','modal');
    return
end
% create path for results
if(exist('plots\Sets\','dir')~=7)
    mkdir('plots\Sets\')
end

%% Parameters
noiseMode={'white','pink'}; % white or pink allowed
numberIterations = 25;
numberClasses = 4;
SNRvalues = [10 20 30 40 50 60 70 80 90 100];
SNRstep = 10;

% preallocate matrices
dataset = zeros(numberIterations,numberClasses,length(SNRvalues),size(algorithms,1));
dataset_ClassMean = zeros(numberIterations,length(SNRvalues),size(algorithms,1));
dataset_ClassIterationMean = zeros(size(algorithms,1),length(SNRvalues));

%% Do plotting
for actualNoiseMode=1:size(noiseMode,2)
    load(['results\evaluation_' noiseMode{actualNoiseMode} '.mat'],'results')
    %% save all NRMSE values in a matrix
    for iteration=1:numberIterations
        iterationName = ['iteration' num2str(iteration)];
        for class=1:numberClasses
            className = ['class' num2str(class)];
            for actualAlgorithm = 1:size(algorithms,1)
                for snrdb = 1:length(SNRvalues)
                    dataset(iteration,class,snrdb,actualAlgorithm)=results.(iterationName).(className)(SNRvalues(snrdb)).(algorithms{actualAlgorithm});
                end
            end
        end
    end
    
    %% calculate mean over classes for each algorithm (for comparison of algorithms over SNR individually)
    for iteration=1:numberIterations
        for snrdb = 1:length(SNRvalues)
            for actualAlgorithm = 1:size(algorithms,1)
                dataset_ClassMean(iteration,snrdb,actualAlgorithm) = mean(dataset(iteration,:,snrdb,actualAlgorithm),'omitnan');
            end
        end
    end
    
    %% calculate mean over classes and iterations for each algorithm (for comparison of algorithms over SNR against each other)
    for snrdb = 1:length(SNRvalues)
        for actualAlgorithm = 1:size(algorithms,1)
            dataset_ClassIterationMean(actualAlgorithm,snrdb) = mean(dataset_ClassMean(:,snrdb,actualAlgorithm),'omitnan');
        end
    end
    
    %% create plots
    SNRlabel = cell(1,length(SNRvalues));
    for snrdb = 1:length(SNRvalues)
        SNRlabel{1,snrdb} = num2str(SNRvalues(snrdb));
    end
    
    for i = 1:size(algorithms,1)
        position(i,:) = (0.8+i/10):1:(10+i/10);
    end
    
    for actualAlgorithm = 1:size(algorithms,1)
            box_O = boxplot(squeeze(dataset_ClassMean(:,:,actualAlgorithm)),'colors',rgb(colors{actualAlgorithm}),'positions',position(actualAlgorithm,:),'width',0.1,'Labels',SNRlabel);
            hold on
            ylabel('NRMSE')
            xlabel('SNR / dB')
            axis([0.9 10.5 0 1])
            out_O = box_O(end,~isnan(box_O(end,:)));
            delete(out_O)
    end
    ax = gca;
    if(mod(size(algorithms,1),2)==0)
        ax.XTick = (position(1,:)+position(end,:))/2;
    else
        ax.XTick = position((1+size(algorithms,1))/2,:);
    end
    c = get(gca, 'Children');
    
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'PaperUnits', 'centimeters');
    currentPos=get(gcf,'Position');
    set(gcf, 'Position', [currentPos(1) currentPos(1) 18 5]);
    
    set(findobj(gcf,'type','axes'),...
        'FontSize', 9,...
        'FontName', 'Times New Roman',...
        'FontWeight','normal',...
        'TitleFontWeight','normal');
    
    a=gcf;
    a.PaperPosition=[0 0 a.Position(3:4)];
    a.PaperSize=a.Position(3:4);
    
    savefig(['plots\Sets\results_' noiseMode{actualNoiseMode} '_' Name 'Boxplot.fig'])
    saveas(gcf,['plots\Sets\results_' noiseMode{actualNoiseMode} '_' Name 'Boxplot.png'])
    matlab2tikz(['plots\Sets\results_' noiseMode{actualNoiseMode} '_' Name 'Boxplot.tex'])
    
    close
    
end
