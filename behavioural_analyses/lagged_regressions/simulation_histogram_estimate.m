clear all;
close all;
clc;

t_stat_fmri  = -3.65;
t_stat_behaviour  = -2.576;

estimate_fmri  = -0.409;
estimate_behaviour  = -0.097;

%X = readtable('data/tstat_permutation.csv');
X = readtable('data/estimate_permutation.csv');

X = table2array(X(:, 1));

%quantiles_tStat = quantile(X,[0.025 0.5 0.975]);
quantiles_estimate = quantile(X,[0.025 0.5 0.975]);

%plot
histogram(X, 20, 'FaceColor',[.75,.75,.75]);
hold on
% get handle to current axes
a = gca;
a.FontSize = 12;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
% set original axes as active
axes(a)
% link axes in case of zooming
linkaxes([a b])
hold on
line([quantiles_estimate(1) quantiles_estimate(1)], [0 20],'color', 'blue');
hold on
line([quantiles_estimate(3) quantiles_estimate(3)], [0 20],'color', 'blue');
hold on
scatter(estimate_behaviour,1,'red'); %insert the result of the lagged analysis from the actually collected data
hold on
scatter(quantiles_estimate(2), 1, 'blue');
%xlabel('t-value for interaction of block type with difference score','fontsize',14);
ylabel('number of group simulations', 'fontsize',14);
hold on
ylim([0,50]);
hold on
xlim([-0.5,0.5]);
hold off

% Compute p
nless = sum(X < estimate_behaviour);
nequal = sum(X == estimate_behaviour);
centile = 100 * (nless + 0.5*nequal) / length(X);
centile = (centile*2)/100

