set(groot, 'defaultaxeslinewidth', 1.25,...
    'defaultlinelinewidth', 1.25, ...
    'defaultstemlinewidth', 1.25, ...
    'defaultstemmarkersize', 2.5, ...,
    'defaultstemmarkerfacecolor', 'auto', ...
    'defaultaxesygrid', 'off', ...
    'defaultaxesyminorgrid', 'off', ...
    'defaultaxesfontsize', 16, ...
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
close all
clear all

fp = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures';
basepath = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode';
save = 0;
%% Level 0.1 speckle

% Ddist analyse filer
% Lfrac = 1/3 her, så omtrent 21 elementer per subarray i full lengde. K_in_lambda = 2
% Resten av verdiene er klassiske. 
fig_ddist_Amp_01 = openfig(fullfile('SpeckleLevel0.1', 'Ddist', 'resolution_AmpCapon.fig'));
fig_ddist_Pow_01 = openfig(fullfile('SpeckleLevel0.1', 'Ddist', 'resolution_PowCapon.fig'));
fig_ddist_DAS_01 = openfig(fullfile('SpeckleLevel0.1', 'Ddist', 'resolution_DAS.fig'));

tmpax = gca(fig_ddist_Pow_01);
tmpax.YLabel.String = 'Amplitude';

%
dist = [0, 0.5, 1, 1.5, 2, 5];

fig_a_EV_dist_01 = figure(100);
fig_a_EV_dist_01.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname      = fullfile('SpeckleLevel0.1', 'ddistEV', ['EV_'  ,num2str(d) ,   'dist_azim_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)

%     clim([-60 0])
    xlim([25 34])
    ylim([1 5])
%     cbb = colorbar()
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d)];
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    if i~=1&&i~=4
        set(ax,'YTickLabel',[]);
    end    
    if i<=3
        set(ax,'XTickLabel',[]);
    end
    
    ax.Parent = tcl;
    ax.Layout.Tile = i;
    set(ax, 'CLim', [-60 5])

    
    i = i+1;
    close(fig)

end
fig_a_EV_dist_01;
sgtitle("Eigenvalues of R through axial line at 0 degrees")

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Label.String = "Eigenvalue [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'Depth [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 14;

fig_d_EV_dist_01 = figure();
fig_d_EV_dist_01.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname      = fullfile('SpeckleLevel0.1', 'ddistEV', ['EV_'  ,num2str(d) ,   'dist_depth_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)

    ylim([1 5])
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d)];
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    if i~=1&&i~=4
        set(ax,'YTickLabel',[]);
    end    
    if i<=3
        set(ax,'XTickLabel',[]);
    end
    
    ax.Parent = tcl;
    ax.Layout.Tile = i;
    set(ax, 'CLim', [-60 5])

    i = i+1;
    close(fig)

end
fig_d_EV_dist_01;
sgtitle("Eigenvalues of R at lateral line through 30[mm]")

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Label.String = "Eigenvalue [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 14;

%% Level 0.25 speckle

% Ddist analyse filer
% Lfrac = 1/3 her, så omtrent 21 elementer per subarray i full lengde. K_in_lambda = 2
% Resten av verdiene er klassiske. 
fig_ddist_Amp_025 = openfig(fullfile('SpeckleLevel0.25', 'Ddist', 'resolution_AmpCapon.fig'));
fig_ddist_Pow_025 = openfig(fullfile('SpeckleLevel0.25', 'Ddist', 'resolution_PowCapon.fig'));
fig_ddist_DAS_025 = openfig(fullfile('SpeckleLevel0.25', 'Ddist', 'resolution_DAS.fig'));

tmpax = gca(fig_ddist_Pow_025);
tmpax.YLabel.String = 'Amplitude';

%
dist = [0, 0.5, 1, 1.5, 2, 5];

fig_a_EV_dist_025 = figure();
fig_a_EV_dist_025.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname      = fullfile('SpeckleLevel0.25', 'ddistEV', ['EV_'  ,num2str(d) ,   'dist_azim_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)

%     clim([-60 0])
    xlim([25 34])
    ylim([1 5])
%     cbb = colorbar()
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d)];
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    if i~=1&&i~=4
        set(ax,'YTickLabel',[]);
    end    
    if i<=3
        set(ax,'XTickLabel',[]);
    end
    
    ax.Parent = tcl;
    ax.Layout.Tile = i;
    set(ax, 'CLim', [-60 5])

    
    i = i+1;
    close(fig)

end
fig_a_EV_dist_025;
sgtitle("Eigenvalues of R through axial line at 0 degrees")

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Label.String = "Eigenvalue [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'Depth [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 14;

fig_d_EV_dist_025 = figure();
fig_d_EV_dist_025.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname      = fullfile('SpeckleLevel0.25', 'ddistEV', ['EV_'  ,num2str(d) ,   'dist_depth_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)

    ylim([1 5])
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d)];
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    if i~=1&&i~=4
        set(ax,'YTickLabel',[]);
    end    
    if i<=3
        set(ax,'XTickLabel',[]);
    end
    
    ax.Parent = tcl;
    ax.Layout.Tile = i;
    set(ax, 'CLim', [-60 5])

    i = i+1;
    close(fig)

end
fig_d_EV_dist_025;
sgtitle("Eigenvalues of R at lateral line through 30[mm]")

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Label.String = "Eigenvalue [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 14;

%% Level 0.5 speckle

% Ddist analyse filer
% Lfrac = 1/3 her, så omtrent 21 elementer per subarray i full lengde. K_in_lambda = 2
% Resten av verdiene er klassiske. 
fig_ddist_Amp_05 = openfig(fullfile('SpeckleLevel0.5', 'Ddist', 'resolution_AmpCapon.fig'));
fig_ddist_Pow_05 = openfig(fullfile('SpeckleLevel0.5', 'Ddist', 'resolution_PowCapon.fig'));
fig_ddist_DAS_05 = openfig(fullfile('SpeckleLevel0.5', 'Ddist', 'resolution_DAS.fig'));

tmpax = gca(fig_ddist_Pow_05);
tmpax.YLabel.String = 'Amplitude';

%
dist = [0, 0.5, 1, 1.5, 2, 5];

fig_a_EV_dist_05 = figure();
fig_a_EV_dist_05.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname      = fullfile('SpeckleLevel0.5', 'ddistEV', ['EV_'  ,num2str(d) ,   'dist_azim_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)

%     clim([-60 0])
    xlim([25 34])
    ylim([1 5])
%     cbb = colorbar()
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d)];
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    if i~=1&&i~=4
        set(ax,'YTickLabel',[]);
    end    
    if i<=3
        set(ax,'XTickLabel',[]);
    end
    
    ax.Parent = tcl;
    ax.Layout.Tile = i;
    set(ax, 'CLim', [-60 5])

    
    i = i+1;
    close(fig)

end
fig_a_EV_dist_05;
sgtitle("Eigenvalues of R through axial line at 0 degrees")

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Label.String = "Eigenvalue [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'Depth [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 14;

fig_d_EV_dist_05 = figure();
fig_d_EV_dist_05.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname      = fullfile('SpeckleLevel0.5', 'ddistEV', ['EV_'  ,num2str(d) ,   'dist_depth_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)

    ylim([1 5])
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d)];
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    if i~=1&&i~=4
        set(ax,'YTickLabel',[]);
    end    
    if i<=3
        set(ax,'XTickLabel',[]);
    end
    
    ax.Parent = tcl;
    ax.Layout.Tile = i;
    set(ax, 'CLim', [-60 5])

    i = i+1;
    close(fig)

end
fig_d_EV_dist_05;
sgtitle("Eigenvalues of R at lateral line through 30[mm]")

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Label.String = "Eigenvalue [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 14;

%% Speckle-level-analysis
d = 2;
speckle_levels = [0.1, 0.25, 0.5];
fig_ev_axial = figure(100);
fig_ev_axial.Position = [100 200 1000 400];
tcl=tiledlayout(1, 4, 'TileSpacing','tight','Padding','compact');
%     genname_az      = "..\Figures\29-Mar-2023\Eigenvalues\Azim\EV_"  +num2str(d)  +   "dist_azim.fig";

genname      = fullfile('29-Mar-2023\Eigenvalues\Azim', ['EV_'  ,num2str(d) ,   'dist_azim.fig']);
fig=openfig(genname);
pause(0.5)
xlim([25 34])
ylim([1 7])
colorbar off
ax = gca;
ax.Title.String = 'Speckle level = 0';
ax.Title.FontSize = 14;
ax.XLabel.Visible = 'off';
ax.YLabel.Visible = 'off';
ax.Parent = tcl;
ax.Layout.Tile = 1;
set(ax, 'CLim', [-60 5])
close(fig)

for i = 1:length(speckle_levels)
    genname      = fullfile(['SpeckleLevel', num2str(speckle_levels(i))], 'ddistEV', ['EV_'  ,num2str(d) ,   'dist_azim_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)
    
    %     clim([-60 0])
    xlim([25 34])
    ylim([1 7])
    %     cbb = colorbar()
    colorbar off
    ax = gca;
    ax.Title.String = ['Speckle level = ', num2str(speckle_levels(i))];
    ax.Title.FontSize = 14;
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    set(ax,'YTickLabel',[]);
    
    ax.Parent = tcl;
    ax.Layout.Tile = i+1;
    set(ax, 'CLim', [-60 5])

    close(fig)
end 
cb = colorbar();
cb.Layout.Tile = 'east';
cb.Title.String = "[dB]";
% cb.Label.Rotation = 270;
% cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])


sgtitle("Eigenvalues of R at axial line through 0 degrees", 'FontSize', 18)

tcl.XLabel.String = 'z [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 16;

%% Lateral line
d = 2;
speckle_levels = [0.1, 0.25, 0.5];
fig_ev_lateral = figure(101);
fig_ev_lateral.Position = [100 200 1000 400];
tcl=tiledlayout(1, 4, 'TileSpacing','tight','Padding','compact');

genname      = fullfile('29-Mar-2023\Eigenvalues\Depth', ['EV_'  ,num2str(d) ,   'dist_depth.fig']);
fig=openfig(genname);
pause(0.5)
ylim([1 7])
colorbar off
ax = gca;
ax.Title.String = 'Speckle level = 0';
ax.Title.FontSize = 14;
ax.XLabel.Visible = 'off';
ax.YLabel.Visible = 'off';
ax.Parent = tcl;
ax.Layout.Tile = 1;
set(ax, 'CLim', [-60 5])
close(fig)

for i = 1:length(speckle_levels)
    genname      = fullfile(['SpeckleLevel', num2str(speckle_levels(i))], 'ddistEV', ['EV_'  ,num2str(d) ,   'dist_depth_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)
    
    ylim([1 7])
    colorbar off
    ax = gca;
    ax.Title.String = ['Speckle level = ', num2str(speckle_levels(i))];
    ax.Title.FontSize = 14;
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    set(ax,'YTickLabel',[]);
    
    ax.Parent = tcl;
    ax.Layout.Tile = i+1;
    set(ax, 'CLim', [-60 5])

    close(fig)
end 
cb = colorbar();
cb.Layout.Tile = 'east';
cb.Title.String = "[dB]";
% cb.Label.Rotation = 270;
% cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])


sgtitle("Eigenvalues of R at lateral line through 30 mm", 'FontSize', 18)

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 16;
%% Resolution/Separability
speckle_levels = [0, 0.1, 0.25, 0.5];
ddist_mat = load(['Data\Fullrun\SpeckleLevel', num2str(speckle_levels(1)) ,'\ddist\ddist_nospeckle.mat']);

tmpmax = ddist_mat.maxDas;

fig_DAS_specklelevels = figure(110);
fig_DAS_specklelevels.Position = [100 200 1000 400];
plot(ddist_mat.x_axis_dmm, db(abs(ddist_mat.DistResDAS(5,:)./tmpmax)), 'DisplayName',['level = ', num2str(speckle_levels(1))])
hold on
for i = 2:length(speckle_levels)
    ddist_mat = load(['Data\Fullrun\SpeckleLevel', num2str(speckle_levels(i)) ,'\ddist\ddist_speckle.mat']);
    plot(ddist_mat.x_axis_dmm, db(abs(ddist_mat.DistResDAS(5,:)./tmpmax)), 'DisplayName',['level = ', num2str(speckle_levels(i))])
    hold on
end
hold off
xlabel('x [mm]')
xlim([-15 15])
ylim([-60 0])
yticks(-60:6:6)
ylabel('Amplitude [dB]')
legend()
title('Separability applying DAS, d=2mm')
grid on

% genname      = fullfile('Figures', '29-Mar-2023', 'Ddist', 'resolution_AmpCapon.fig');
ddist_mat = load(['Data\Fullrun\SpeckleLevel', num2str(speckle_levels(1)) ,'\ddist\ddist_nospeckle.mat']);

fig_ampCap_specklelevels = figure(111);
fig_ampCap_specklelevels.Position = [100 200 1000 400];
plot(ddist_mat.x_axis_dmm, db(abs(ddist_mat.DistResMVAmp(5,:)./tmpmax)), 'DisplayName',['level = ', num2str(speckle_levels(1))])
hold on
for i = 2:length(speckle_levels)
    ddist_mat = load(['Data\Fullrun\SpeckleLevel', num2str(speckle_levels(i)) ,'\ddist\ddist_speckle.mat']);
    plot(ddist_mat.x_axis_dmm, db(abs(ddist_mat.DistResMVAmp(5,:)./tmpmax)), 'DisplayName',['level = ', num2str(speckle_levels(i))])
    hold on
end
hold off
xlabel('x [mm]')
xlim([-15 15])
ylim([-60 0])
yticks(-60:6:6)
ylabel('Amplitude [dB]')
legend()
title('Separability applying Amplitude Capon, d=2mm')
grid on

ddist_mat = load(['Data\Fullrun\SpeckleLevel', num2str(speckle_levels(1)) ,'\ddist\ddist_nospeckle.mat']);

fig_powCap_specklelevels = figure(112);
fig_powCap_specklelevels.Position = [100 200 1000 400];
plot(ddist_mat.x_axis_dmm, 10*log10(abs(ddist_mat.DistResMVPow(5,:))), 'DisplayName',['level = ', num2str(speckle_levels(1))])
hold on
for i = 2:length(speckle_levels)
    ddist_mat = load(['Data\Fullrun\SpeckleLevel', num2str(speckle_levels(i)) ,'\ddist\ddist_speckle.mat']);
    plot(ddist_mat.x_axis_dmm, 10*log10(abs(ddist_mat.DistResMVPow(5,:))), 'DisplayName',['level = ', num2str(speckle_levels(i))])
    hold on
end
hold off
xlabel('x [mm]')
xlim([-15 15])
ylim([-60 -6])
yticks(-60:6:6)
ylabel('Amplitude [dB]')
legend()
title('Separability applying Power Capon, d=2mm')
grid on


%%
if save
    exportgraphics(fig_DAS_specklelevels, fullfile(basepath, 'Figures', 'Overleaf', 'SpeckleLevelAnalysis', 'DAS_speckleLevel_resolution.pdf'))
    exportgraphics(fig_ampCap_specklelevels, fullfile(basepath, 'Figures', 'Overleaf', 'SpeckleLevelAnalysis', 'AmpCap_speckleLevel_resolution.pdf'))
    exportgraphics(fig_powCap_specklelevels, fullfile(basepath, 'Figures', 'Overleaf', 'SpeckleLevelAnalysis', 'PowCap_speckleLevel_resolution.pdf'))
    exportgraphics(fig_ev_lateral, fullfile(basepath, 'Figures', 'Overleaf', 'SpeckleLevelAnalysis', 'EV_speckleLevel_lateral.pdf'))
    exportgraphics(fig_ev_axial, fullfile(basepath, 'Figures', 'Overleaf', 'SpeckleLevelAnalysis', 'EV_speckleLevel_axial.pdf'))
end



