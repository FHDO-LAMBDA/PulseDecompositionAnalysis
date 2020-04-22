clc
clear all
close all

% NaN is omitted


%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% load file containing Algorithmnames here
load('algorithmsCOMPLETE2.mat','algorithms');
% load file specifying colors for the plots here
load('colorsCOMPLETE2.mat','colors');
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
% create path for plots
if(exist('plots\','dir')~=7)
    mkdir('plots\')
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
dataset_mean = zeros(1,size(algorithms,1));
dataset_median = zeros(1,size(algorithms,1));

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
    
    %% calculate and save mean NRMSE in table
    % mean NRMSE over iterations and SNR and classes for each algorithm
    for actualAlgorithm = 1:size(algorithms,1)
        dataset_mean(1,actualAlgorithm) = mean(dataset_ClassIterationMean(actualAlgorithm,:),'omitnan');
        dataset_median(1,actualAlgorithm) = median(dataset_ClassIterationMean(actualAlgorithm,:),'omitnan');
    end
    % store results in table
    dataset_meanTable = array2table(dataset_mean,'VariableNames',algorithms);
    dataset_medianTable = array2table(dataset_median,'VariableNames',algorithms);
    save(['results\NRMSEtable' noiseMode{actualNoiseMode} '.mat'],'dataset_meanTable','dataset_medianTable');
    
    %% create plots
    % compare algorithms against one another
    fig0 = figure('Name','CompareAllAlgorithms');
    hold on
    for actualAlgorithm = 1:size(algorithms,1)
        plot(min(SNRvalues):SNRstep:max(SNRvalues),...
            dataset_ClassIterationMean(actualAlgorithm,:),...
            'Color',rgb(colors{actualAlgorithm}))
    end
    ylabel('NRMSE_{deriv2}')
    xlabel('SNR/dB')
    legend(algorithms','Location','eastoutside')
    axis([min(SNRvalues) max(SNRvalues) -10 1])
    savefig(['plots\CompareAllAlgorithms_' noiseMode{actualNoiseMode} '.fig'])
    saveas(fig0,['plots\CompareAllAlgorithms_' noiseMode{actualNoiseMode} '.png'])
    matlab2tikz(['plots\CompareAllAlgorithms_' noiseMode{actualNoiseMode} '.tex'])
    
    % compare algorithms individually
    SNRlabel = cell(1,length(SNRvalues));
    for snrdb = 1:length(SNRvalues)
        SNRlabel{1,snrdb} = num2str(SNRvalues(snrdb));
    end
    
    for actualAlgorithm = 1:size(algorithms,1)
        currentMatrix = squeeze(dataset_ClassMean(:,:,actualAlgorithm));
        figure('Name',algorithms{actualAlgorithm});
        boxplot(currentMatrix,'Labels',SNRlabel)
        ylabel('NRMSE_{deriv2}')
        xlabel('SNR/dB')
        axis([-inf inf 0 1])
        if(exist(['plots\' algorithms{actualAlgorithm} '\'],'dir')~=7)
            mkdir(['plots\' algorithms{actualAlgorithm} '\'])
        end
        savefig(['plots\' algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '.fig'])
        saveas(gcf,['plots\' algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '.png'])
        matlab2tikz(['plots\' algorithms{actualAlgorithm} '\' noiseMode{actualNoiseMode} '.tex'])
    end
        
end