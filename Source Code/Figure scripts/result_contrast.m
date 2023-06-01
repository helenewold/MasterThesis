close all
clear all

%%
set(groot, 'defaultaxeslinewidth', 1.25,...
    'defaultlinelinewidth', 1.25, ...
    'defaultstemlinewidth', 1.25, ...
    'defaultstemmarkersize', 2.5, ...,
    'defaultstemmarkerfacecolor', 'auto', ...
    'defaultaxesygrid', 'off', ...
    'defaultaxesyminorgrid', 'off', ...
    'defaultfigureposition', [50 50 600 600], ...
    'defaultaxesfontsize', 18, ...
    'defaultaxesyminortick', 'on', ...
    'defaultaxestickdir', 'out', ...
    'defaultAxesTickDirMode', 'manual', ...
    'defaultaxesticklength', 6*[0.001 0.001], ...
    'defaulttextinterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ... 
    'defaultfigurecolor', 'white' ...
    );
%%
save = 0

data = load("Data\Contrast\contrast_metrics.mat")
%%
CNR_DAS_circ = data.CNR_DAS_circ;
CR_DAS_circ = data.CR_DAS_circ;
gCNR_DAS_circ = data.gCNR_DAS_circ;

CNR_DAS_orig = data.CNR_DAS_orig;
CR_DAS_orig = data.CR_DAS_orig;
gCNR_DAS_orig = data.gCNR_DAS_orig;

CNR_DAS_small = data.CNR_DAS_small;
CR_DAS_small = data.CR_DAS_small;
gCNR_DAS_small = data.gCNR_DAS_small;

%%
CNR_circ = data.CNR_circ;
CR_circ = data.CR_circ;
gCNR_circ = data.gCNR_circ;

CNR_orig = data.CNR_orig;
CR_orig = data.CR_orig;
gCNR_orig = data.gCNR_orig;

CNR_small = data.CNR_small;
CR_small = data.CR_small;
gCNR_small = data.gCNR_small;

Lelms = data.Lelms;
DL = data.DL;
%%
partLelm = round([1/4, 1/3, 1/2]*64);
fig1_CNR = figure(1);
% fig1_CNR.Position = [100 200 750 450];

tcl=tiledlayout(3, 1, 'TileSpacing','tight','Padding','compact');

tmpfig = figure();

plot(DL*100, CNR_circ(:, partLelm))
hold on
yline(CNR_DAS_circ, '--r')
hold off
ylabel("CNR")
title("r = 5 mm")
xlabel("Diagonal Load \%")
% ylim([0.838 CNR_DAS_small])

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 3;
% set(ax, 'XTickLabel',[]);

plot(DL*100, CNR_orig(:, partLelm))
hold on
yline(CNR_DAS_orig, '--r')
hold off
ylabel("CNR")
title("r = 3.5 mm")
% ylim([0.838 CNR_DAS_small])

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 2;
set(ax, 'XTickLabel',[]);

plot(DL*100, CNR_small(:, partLelm))
hold on
yline(CNR_DAS_small, '--r')
hold off
ylabel("CNR")
title("r = 2.5 mm")
legend('L = 1/4', 'L = 1/3', 'L = 1/2', 'DAS', 'Location','southeast')
legend('Orientation', 'horizontal', 'Location', 'northoutside')
% ylim([0.838 CNR_DAS_small])

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 1;
set(ax, 'XTickLabel',[]);

close(tmpfig)
%% gCNR
fig1_CR = figure(2);
% fig1_CR.Position = [100 200 750 450];

tcl=tiledlayout(3, 1, 'TileSpacing','tight','Padding','compact');

tmpfig = figure();

plot(DL*100, CR_circ(:, partLelm))
hold on
yline(CR_DAS_circ, '--r')
hold off
title("r = 5 mm")
ylabel("CR")
xlabel("Diagonal Load \%")

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 3;
ax.YDir = 'reverse';

plot(DL*100, CR_orig(:, partLelm))
hold on
yline(CR_DAS_orig, '--r')
hold off
title("r = 3.5 mm")
ylabel("CR")

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 2;
set(ax, 'XTickLabel',[]);
ax.YDir = 'reverse';

plot(DL*100, CR_small(:, partLelm))
hold on
yline(CR_DAS_small, '--r')
hold off
title("r = 2.5 mm")
ylabel("CR")
legend('L = 1/4', 'L = 1/3', 'L = 1/2', 'DAS', 'Location','southeast')
legend('Orientation', 'horizontal', 'Location', 'northoutside')
ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 1;
set(ax, 'XTickLabel',[]);
ax.YDir = 'reverse';

close(tmpfig)
%%
fig1_gCNR = figure(3);
% fig1_gCNR.Position = [100 200 750 450];

tcl=tiledlayout(3, 1, 'TileSpacing','tight','Padding','compact');

tmpfig = figure();
plot(DL*100, gCNR_circ(:, partLelm))
hold on
yline(gCNR_DAS_circ, '--r')
hold off
title("r = 5 mm")

xlabel("Diagonal Load \%")
ylabel("gCNR")

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 3;
% set(ax, 'XTickLabel',[]);

plot(DL*100, gCNR_orig(:, partLelm))
hold on
yline(gCNR_DAS_orig, '--r')
hold off
title('r = 3.5 mm')
ylabel("gCNR")

% xlabel("Diagonal Load %")

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 2;
set(ax, 'XTickLabel',[]);

plot(DL*100, gCNR_small(:, partLelm))
hold on
yline(gCNR_DAS_small, '--r')
hold off
title("r = 2.5 mm")

% xlabel("Diagonal Load %")
ylabel("gCNR")
legend('L = 1/4', 'L = 1/3', 'L = 1/2', 'DAS', 'Location','southeast')
legend('Orientation', 'horizontal', 'Location', 'northoutside')
ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 1;
set(ax, 'XTickLabel',[]);


close(tmpfig)
%%
DLs = [1 2 8 10 15 33];
fig2_CNR = figure(4);
% fig2_CNR.Position = [100 200 750 450];


tcl=tiledlayout(3, 1, 'TileSpacing','tight','Padding','compact');

tmpfig = figure();


plot(Lelms, CNR_circ(DLs, 2:end))
hold on
yline(CNR_DAS_circ, '--r')
hold off
ylabel("CNR")
xlabel("Subarray size")
xlim([2 63])
% ylim([0.6 0.9])
title("r = 5 mm")

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 3;
% set(ax, 'XTickLabel',[]);


plot(Lelms, CNR_orig(DLs, 2:end))
hold on
yline(CNR_DAS_orig, '--r')
hold off
title("r = 3.5 mm")
ylabel("CNR")
xlim([2 63])

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 2;
set(ax, 'XTickLabel',[]);
% ax.YDir = 'reverse';

plot(Lelms, CNR_small(DLs, 2:end))
hold on
yline(CNR_DAS_small, '--r')
hold off
title("r = 2.5 mm")

ylabel("CNR")
xlim([2 63])
ylim([0.6 1.001])

legend('DL=0\%', 'DL=0.0001\%', 'DL=0.1\%','DL=1\%','DL=10\%', 'DL=100\%', 'DAS')
legend('NumColumns', 4, 'Orientation', 'horizontal', 'Location', 'northoutside')

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 1;
set(ax, 'XTickLabel',[]);

close(tmpfig)


%%
fig2_CR = figure(5);
% fig2_CR.Position = [100 200 750 450];


tcl=tiledlayout(3, 1, 'TileSpacing','tight','Padding','compact');

tmpfig = figure();

plot(Lelms, CR_circ(DLs, 2:end))
hold on
yline(CR_DAS_circ, '--r')
hold off
ylabel("CR")
xlim([2 63])
title("r = 5mm")
% ylim([0.6 0.9])

xlabel("Subarray size")
ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 3;
% set(ax, 'XTickLabel',[]);
ax.YDir = 'reverse';


plot(Lelms, CR_orig(DLs, 2:end))
hold on
yline(CR_DAS_orig, '--r')
hold off
ylabel("CR")
title("r = 3.5 mm")
xlim([2 63])

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 2;
set(ax, 'XTickLabel',[]);
ax.YDir = 'reverse';


plot(Lelms, CR_small(DLs, 2:end))
hold on
yline(CR_DAS_small, '--r')
hold off
ylabel("CR")
xlim([2 63])
title("r = 2.5mm")

legend('DL=0\%', 'DL=0.0001\%', 'DL=0.1\%','DL=1\%','DL=10\%', 'DL=100\%', 'DAS')
legend('NumColumns', 4, 'Orientation', 'horizontal', 'Location', 'northoutside')

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 1;
set(ax, 'XTickLabel',[]);
ax.YDir = 'reverse';

close(tmpfig)

%%
fig2_gCNR = figure(6);
% fig2_gCNR.Position = [100 200 750 450];


tcl=tiledlayout(3, 1, 'TileSpacing','tight','Padding','compact');

tmpfig = figure();

plot(Lelms, gCNR_circ(DLs, 2:end))
hold on
yline(gCNR_DAS_circ, '--r')
hold off
ylabel("gCNR")
xlabel("Subarray size")
title("r = 5mm")
xlim([2 63])
ylim([0.73 0.85])

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 3;


plot(Lelms, gCNR_orig(DLs, 2:end))
hold on
yline(gCNR_DAS_orig, '--r')
hold off
title("r = 3.5mm")
ylabel("gCNR")
ylim([0.96 1.001])

xlim([2 63])


% xlabel("Diagonal Load %")


ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 2;
set(ax, 'XTickLabel',[]);


plot(Lelms, gCNR_small(DLs, 2:end))
hold on
yline(gCNR_DAS_small, '--r')
hold off
% xlabel("Subarray size")
ylabel("gCNR")
title("r = 2.5mm")

ylim([0.99 1.001])
xlim([2 63])

legend('DL=0\%', 'DL=0.0001\%', 'DL=0.1\%','DL=1\%','DL=10\%', 'DL=100\%', 'DAS')
legend('NumColumns', 4, 'Orientation', 'horizontal', 'Location', 'northoutside')

ax = gca;
ax.Parent = tcl;
ax.Layout.Tile = 1;
set(ax, 'XTickLabel',[]);


close(tmpfig)



%%
if save
    path = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\Overleaf\Contrast';
%     mkdir(path);
    exportgraphics(fig1_CR, fullfile(path, 'CR_x_axis_DL.pdf'))
    exportgraphics(fig2_CR, fullfile(path, 'CR_x_axis_Lelm.pdf'))    
    exportgraphics(fig1_CNR, fullfile(path, 'CNR_x_axis_DL.pdf'))
    exportgraphics(fig2_CNR, fullfile(path, 'CNR_x_axis_Lelm.pdf'))    
    exportgraphics(fig1_gCNR, fullfile(path, 'gCNR_x_axis_DL.pdf'))
    exportgraphics(fig2_gCNR, fullfile(path, 'gCNR_x_axis_Lelm.pdf'))
end