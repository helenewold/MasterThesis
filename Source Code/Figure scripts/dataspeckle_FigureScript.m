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
save = 0;
cd("C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Source Code")
basepath = "C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode";
%% D dist analysis
ddist_analysis = open(fullfile('Data', 'Fullrun', 'SpeckleLevel0.25', 'ddist', 'ddist_speckle.mat'));
distDAS     = ddist_analysis.DistResDAS;
distAMP     = ddist_analysis.DistResMVAmp;
distPOW     = ddist_analysis.DistResMVPow;
azimuth_ax  = ddist_analysis.azimuth_MLA;
depth_ax    = ddist_analysis.depth_axis;
maxDas      = ddist_analysis.maxDasVals;
x_axis_dmm  = ddist_analysis.x_axis_dmm;
%% Plotting
dist = [0, 0.5, 1, 1.5, 2, 5];

L = {['d=',num2str(dist(1)), 'mm'],...
    ['d=',num2str(dist(2)), 'mm'],...
    ['d=',num2str(dist(3)), 'mm'],...
    ['d=',num2str(dist(4)), 'mm'],...
    ['d=',num2str(dist(5)), 'mm'],...
    ['d=',num2str(dist(6)), 'mm']
    };
lt = {'-', ':', ':', '--', '--','-'};

figdistDas = figure(1);
drawnow;
figdistDas.Position = [100 200 1000 400];
for iii = 1:length(dist)
    drawnow
    plot(x_axis_dmm, db(abs(distDAS(iii,:)./maxDas(end))), 'LineWidth',1.25, 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii)))
    hold on
end
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm] - DAS")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 5])
xlim([-5 10])
yticks(-60:6:6)
grid on
% %%
figdistMVAmp = figure(2);
drawnow;
figdistMVAmp.Position = [100 200 1000 400];
for iii = 1:length(dist)
    drawnow
    plot(x_axis_dmm, db(abs(distAMP(iii,:)./maxDas(end))), 'LineWidth',1.25, 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii)))
    hold on
end
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm] - Amplitude Capon")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 5])
xlim([-5 10])
yticks(-60:6:6)
grid on

% %%
figdistMVPow = figure(3);
drawnow;
figdistMVPow.Position = [100 200 1000 400];
for iii = 1:length(dist)
    drawnow
    plot(x_axis_dmm, 10*log10(abs(distPOW(iii,:))), 'LineWidth',1.25, 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii)))
    hold on
end
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm] - Power Capon")
xlabel("x [mm]")
ylabel("Amplitude")
ylim([-60 5])
xlim([-5 10])
yticks(-120:6:6)
grid on

if save
    path1 = fullfile(basepath, 'Figures', 'Overleaf');
    exportgraphics(figdistDas, fullfile(path1, 'SpeckleDdist', 'resolution_speckle_distDAS.pdf'))
    exportgraphics(figdistMVAmp, fullfile(path1, 'SpeckleDdist', 'resolution_speckle_distAMP.pdf'))
    exportgraphics(figdistMVPow, fullfile(path1, 'SpeckleDdist', 'resolution_speckle_distPOW.pdf'))
end
%% Eigenvalues fra analysen over
ev_ddist_path = fullfile(basepath, 'Figures', 'Fullrun', 'SpeckleLevel0.25','ddistEV');
fig_EV_ddist_az = figure(4);
fig_EV_ddist_az.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');

for i = 1:6
    genname = fullfile(ev_ddist_path, ['EV_', num2str(dist(i)), 'dist_azim_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)

    ylim([1 5])
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(dist(i))];
    ax.Title.FontSize = 14;
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

%     i = i+1;
    close(fig)
end
sgtitle("Eigenvalues of R through axial line at 0 degrees", 'FontSize', 18)

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Title.String = "[dB]";
% cb.Label.Rotation = 270;
% cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 16;

fig_EV_ddist_dp = figure(5);
fig_EV_ddist_dp.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');

for i = 1:6
    genname = fullfile(ev_ddist_path, ['EV_', num2str(dist(i)), 'dist_depth_speckle.fig']);
    fig=openfig(genname);
    pause(0.5)

    ylim([1 7])
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(dist(i))];
    ax.Title.FontSize = 14;
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

%     i = i+1;
    close(fig)
end
sgtitle("Eigenvalues of R at lateral line through 30[mm]", 'FontSize', 18)

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Title.String = "[dB]";
% cb.Label.Rotation = 270;
% cb.Label.Position(1) = 3;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 16;
if save
    path1 = fullfile(basepath, 'Figures', 'Overleaf');
    exportgraphics(fig_EV_ddist_az, fullfile(path1, 'SpeckleDdist', 'EV_ddist_speckle_az.pdf'))
    exportgraphics(fig_EV_ddist_dp, fullfile(path1, 'SpeckleDdist', 'EV_ddist_speckle_dp.pdf'))
end
%% Subarray averaging analysis (Kjører fullt for øyeblikket......)
L_elm_analysis = open(fullfile(basepath, 'Data', 'Fullrun', 'SpeckleLevel0.25', 'CNDLLelm', ['CN_DL_Lelm_speckle.mat']));

CN_results = L_elm_analysis.CN_results;
DLs = L_elm_analysis.DL;
DasRes = L_elm_analysis.DasRes;
TotRes = L_elm_analysis.TotRes;
TotRes_Pow = L_elm_analysis.TotRes_Pow;
maxDas_Lelm = L_elm_analysis.maxDas;
x_axis_dmm = L_elm_analysis.x_axis_dmm;

x_axis_dmm = linspace(x_axis_dmm(1), x_axis_dmm(end), size(TotRes, 3))

%% Condition number
fig_CNDLLelm = figure(10);
fig_CNDLLelm.Position = [200 100 1000 500];
Ls = {'max', 'midpoint', 'leftpoint', 'rightpoint', 'mid-5', 'mid+5'};
for i = 1:6
    tmpinds = find(CN_results(i,:,:) > 1e16);
    CN_results(i,tmpinds) = 0;
    subplot(2,3,i)
    imagesc(log10(squeeze(CN_results(i,:,:))))
    yticks(1:size(CN_results,2))
    yticklabels(split(num2str(DLs*100)))
    xlabel("Subarray size")
    ylabel("Diagonal load (\%)")
    title(string(Ls(i)))
    c = colorbar();
    c.Label.String = "Condition number";
end

figure_DL_CN = figure(11);
figure_DL_CN.Position = [200 100 750 400];
for ij = 2:10:53
    Ls = ['Lelm=', num2str(ij)];
    plot(DLs*100, log10(CN_results(2,:,ij)), 'DisplayName',string(Ls), 'LineWidth',1)
    hold on
end
hold off
legend('Location', 'northeast')
xlabel('Diagonal Load (\%)')
ylabel('log10(CN)')
title('Condition number at point between two scatterers')

Ls = {'max', 'midpoint', 'leftpoint', 'rightpoint', 'mid-5', 'mid+5'};


figure_Lelm_CN = figure(12);
figure_Lelm_CN.Position = [200 100 400 200];
for ij = 1:6
%     Ls = ['Lelm=', num2str(ij)];
    plot(1:63, log10(squeeze(CN_results(ij,10,:))), 'DisplayName',string(Ls(ij)), 'LineWidth',1)
    hold on
end
hold off
legend('Location', 'southeast')
xlabel('Subarray Size')
ylabel('log10(CN)')
xlim([2 63])
title('Condition number for different locations around scatterers')

if save
    path2 = fullfile(basepath, 'Figures', 'Overleaf', 'Speckle_CNDLLelm');
    exportgraphics(fig_CNDLLelm,    fullfile(path2, 'CN_DL_Lelm_total.pdf'))
    exportgraphics(figure_DL_CN,    fullfile(path2, 'CN_DL_mid.pdf'))
    exportgraphics(figure_Lelm_CN,  fullfile(path2, 'CN_Lelm_points.pdf'))
end

%% Resolution subarray size
% TotRes(4,1,:) = DasRes(:);
TotRes(10,1,:) = DasRes(:);
fig_full_subarr = figure(13);
fig_full_subarr.Position = [200 100 500 500];
imagesc(x_axis_dmm, 1:63, db(abs(squeeze(TotRes(10,:,:)./maxDas_Lelm))))
xlim([-10 10])
cb = colorbar();
cb.Label.String = "Amplitude [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 4;
xlabel('x[mm]')
ylabel('Subarray size')
title('Lateral line through 30[mm]')

TotResNorm = zeros(63, 384);
for i = 1:63
    TotResNorm(i,:) = squeeze(TotRes(10,i,:)./max(TotRes(10,i,:)));
end

fig_full_subarr_norm = figure(14);
fig_full_subarr_norm.Position = [200 100 500 500];
imagesc(x_axis_dmm, 1:63, db(abs(TotResNorm)))
xlim([-10 10])
cb = colorbar();
cb.Label.String = "Amplitude [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 4;
xlabel('x[mm]')
ylabel('Subarray size')
title('Lateral line through 30[mm]')
 %%
 Lelms = 2:10:53;
L = {'DAS',...
    ['Lelm=',num2str(Lelms(1))],...
    ['Lelm=',num2str(Lelms(2))],...
    ['Lelm=',num2str(Lelms(3))],...
    ['Lelm=',num2str(Lelms(4))],...
    ['Lelm=',num2str(Lelms(5))],...
    ['Lelm=',num2str(Lelms(6))]
    };
lt = {'-', ':', ':', '--', '--','-'};

LelmRes = zeros(6,128*3);
LelmRes(1,:) = DasRes(:);
LelmRes(2:6,:) = squeeze(TotRes(10,2:10:51,:));

fig_part_subarr = figure(15);
fig_part_subarr.Position = [200 100 750 400];
plot(x_axis_dmm, db(abs(LelmRes(1,:)./maxDas_Lelm)), '--', 'LineWidth',1.25, 'DisplayName', string(L(1)))
hold on
for iii = 2:6
    drawnow
    plot(x_axis_dmm, db(abs(LelmRes(iii,:)./maxDas_Lelm)), 'LineWidth',1.25, 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii)))
    hold on
end
% drawnow
hold off
title('Lateral line through 30 mm')
legend('Location', 'southeast', 'NumColumns',2)
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 5])
xlim([-5 5])
yticks(-120:6:6)
grid on

%%
TotResPOWNorm = zeros(63, 384);
for i = 1:63
    TotResPOWNorm(i,:) = squeeze(TotRes_Pow(10,i,:)./max(TotRes_Pow(10,i,:)));
end

fig_full_subarr_Pow = figure(16);
fig_full_subarr_Pow.Position = [200 100 1000 500];
imagesc(x_axis_dmm, 2:63, 10*log10(abs(TotResPOWNorm(2:end,:))))
cb = colorbar();
cb.Label.String = "Amplitude";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 4;
xlabel('x[mm]')
ylabel('Subarray size')
title('Lateral line through 30[mm] - Power Capon')
LelmRes_Pow = squeeze(TotRes_Pow(10, 2:10:53,:));

%%
fig_part_subarr_Pow = figure(17);
fig_part_subarr_Pow.Position = [200 100 750 400];
% plot(x_axis_dmm, 10*log10(abs(TotRes_Pow(1,:))), '--','DisplayName',string(L(1)), 'LineWidth',1.25)
hold on
for iii = 1:6
%     drawnow
    plot(x_axis_dmm, 10*log10(abs(LelmRes_Pow(iii,:))), 'LineWidth',1.25, 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii+1)))
    hold on
end
% drawnow
hold off
title('Lateral line through 30 mm - Power Capon')
legend('Location', 'southeast', 'NumColumns',3)
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 5])
xlim([-15 15])
yticks(-120:6:6)
grid on
%%
if save
    path1 = fullfile(basepath, 'Figures', 'Overleaf', 'SpeckleSubarrayAnalysis');
    exportgraphics(fig_full_subarr,         fullfile(path1, 'Lelm_res_speckle_2D.pdf'))
    exportgraphics(fig_full_subarr_norm,    fullfile(path1, 'Lelm_res_speckle_2D_norm.pdf'))
    exportgraphics(fig_part_subarr,         fullfile(path1, 'Lelm_res_speckle_part.pdf'))
    exportgraphics(fig_full_subarr_Pow,     fullfile(path1, 'Pow_Lelm_res_speckle_2D.pdf'))
    exportgraphics(fig_part_subarr_Pow,     fullfile(path1, 'Pow_Lelm_res_speckle_part.pdf'))
    
end
