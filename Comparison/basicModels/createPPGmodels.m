clc 
close all
clear all

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\..\Algorithms');
% load file containing Algorithmnames here
load('algorithmsCOMPLETE2.mat','algorithms');
% check that there is an algorithm for every entry in the list
for actualAlgorithm = 1:size(algorithms,1)
    if(~(isfile(['..\..\Algorithms\' 'decompose' algorithms{actualAlgorithm} '.m'])))
        errordlg('Unknown decomposition algorithm selected','Input Error','modal');
        return
    end
end
% specify source and results folder
sourceFolder = 'Reference\';
resultsFolder = 'Parameter\';

%% Parameters
samplingFrequency = 1000;
beatStart = 1;
beatEnd = 700;

%% Create model PPGs
for class = 1:4
    load(['Reference\Parameter_',num2str(class),'_gauss4.mat'],'signal'); % load reference signal for current class
    signal = signal(beatStart:beatEnd); % cut reference beat from 1:700
    trenddata = interp1([beatStart beatEnd], [signal(beatStart) signal(beatEnd)],1:length(signal),'linear','extrap'); % calculate straight line from first to last point
    signal = signal - trenddata; % do detrending by removing the straight line
    signal = (signal - min(signal))/(max(signal) - min(signal)); % normalize reference beat
    t_ppg = 0:1/samplingFrequency:(length(signal)-1)/samplingFrequency;
    for actualAlgorithm = 1:size(algorithms,1)
        decomposeFuncName = ['decompose' algorithms{actualAlgorithm}];
        [signal_mod,y,opt_values] = feval(decomposeFuncName,signal,samplingFrequency);
        figure('Name',[algorithms{actualAlgorithm,1} ' class ' num2str(class)])
        hold on
        kernels = plot(repmat(t_ppg,size(y,1),1)',cell2mat(y)',':k');
        model = plot(t_ppg',signal_mod','k');
        real = plot(t_ppg',signal','r');
        hold off
        legend([kernels(1); model; real],'basic components','model beat','original beat');
        legend('boxoff');
        xlabel('time/s');
        ylabel('amplitude/a.u.');
        xlim([t_ppg(1) t_ppg(end)]);
        ylim([0 1.1]);
        if(exist([resultsFolder algorithms{actualAlgorithm} '\'],'dir')~=7)
            mkdir([resultsFolder algorithms{actualAlgorithm} '\'])
        end
        save([resultsFolder algorithms{actualAlgorithm} '\' num2str(class) '.mat'],'signal_mod','y','opt_values')
        savefig([resultsFolder algorithms{actualAlgorithm} '\' num2str(class) '.fig'])
    end
end