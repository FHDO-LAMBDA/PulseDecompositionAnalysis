clear all
close all
clc

scalingFactor=1;
myFontType='Times New Roman';
myFontSize=6/scalingFactor;
myFigureSize=[2 4];

load('Parameter_1_gauss4.mat')
signal = signal(1:700); % cut reference beat from 1:700
t = 0:(1/1000):(length(signal)-1)/1000;
plot(t,signal,'k','LineWidth',1)
xticks([0 0.5])
xticklabels({'0','0.5'})
xlabel('time / s','FontSize',myFontSize)
ylabel('amplitude / a.u.','FontSize',myFontSize)
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
savefig(gcf,'Class1.fig')
saveas(gcf,'Class1.png')
matlab2tikz('Class1.tex')

load('Parameter_2_gauss4.mat')
signal = signal(1:700); % cut reference beat from 1:700
t = 0:(1/1000):(length(signal)-1)/1000;
plot(t,signal,'k','LineWidth',1)
xticks([0 0.5])
xticklabels({'0','0.5'})
xlabel('time / s','FontSize',myFontSize)
ylabel('amplitude / a.u.','FontSize',myFontSize)
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
savefig(gcf,'Class2.fig')
saveas(gcf,'Class2.png')
matlab2tikz('Class2.tex')

load('Parameter_3_gauss4.mat')
signal = signal(1:700); % cut reference beat from 1:700
t = 0:(1/1000):(length(signal)-1)/1000;
plot(t,signal,'k','LineWidth',1)
xticks([0 0.5])
xticklabels({'0','0.5'})
xlabel('time / s','FontSize',myFontSize)
ylabel('amplitude / a.u.','FontSize',myFontSize)
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
savefig(gcf,'Class3.fig')
saveas(gcf,'Class3.png')
matlab2tikz('Class3.tex')

load('Parameter_4_gauss4.mat')
signal = signal(1:700); % cut reference beat from 1:700
t = 0:(1/1000):(length(signal)-1)/1000;
plot(t,signal,'k','LineWidth',1)
xticks([0 0.5])
xticklabels({'0','0.5'})
xlabel('time / s','FontSize',myFontSize)
ylabel('amplitude / a.u.','FontSize',myFontSize)
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
savefig(gcf,'Class4.fig')
saveas(gcf,'Class4.png')
matlab2tikz('Class4.tex')