clc
close all
clear all


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
% create path for results
if(exist('plots\','dir')~=7)
    mkdir('plots\')
end

%% Parameters
movArt = {'mov1','mov2'}; % mov1 and mov2 allowed
numberClasses = 4;

% preallocate matrices
dataset = zeros(size(algorithms,1),numberClasses);

%% Do plotting
for actualMovArt = 1:size(movArt,2)
    load(['results\evaluation_' movArt{actualMovArt} '.mat'])
    
    
    for class=1:numberClasses
        className = ['class' num2str(class)];
        for actualAlgorithm = 1:size(algorithms,1)
            dataset(actualAlgorithm,class)=results.(className).(algorithms{actualAlgorithm});
        end
    end
    

    figure
    hold on
    for actualAlgorithm = 1:size(algorithms,1)
        plot([1 2 3 4],dataset(actualAlgorithm,:),'Color',rgb(colors{actualAlgorithm}))
    end
    ylabel('NRMSE_{deriv2}')
    xlabel('class')
    axis([1 4 0 1])
    ax = gca;
    ax.XTick = [1 2 3 4];
    legend(algorithms','Location','eastoutside')
    savefig(['plots\CompareAllAlgorithms_' movArt{actualMovArt} '.fig'])
    saveas(gcf,['plots\CompareAllAlgorithms_' movArt{actualMovArt} '.png'])
    matlab2tikz(['plots\CompareAllAlgorithms_' movArt{actualMovArt} '.tex'])
    
end