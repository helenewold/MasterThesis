% Tanke bak scriptet: Ha lagret alle figurer som fig-fil, for så å bruke
% dette scriptet til å kunne laste inn figuren. Greit dersom det trengs
% bittesmå endringer på et script...
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
%% Figur subarray length

path = fullfile('C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode', 'Data', '02-Apr-2023');
filename = 'Lelm_analysis2.mat';

Lelmanalysis = open(fullfile(path, filename));

LelmRes = Lelmanalysis.LelmRes;
Lelms = Lelmanalysis.Lelms;
TotRes = Lelmanalysis.TotRes;
maxDas = Lelmanalysis.maxDas;
x_axis_dmm = Lelmanalysis.x_axis_dmm;

LelmRes_Norm(:,:) = abs(LelmRes(:,:))./maxDas;
%
L = {'DAS',...
    ['Lelm=',num2str(Lelms(1))],...
    ['Lelm=',num2str(Lelms(2))],...
    ['Lelm=',num2str(Lelms(3))],...
    ['Lelm=',num2str(Lelms(4))],...
    ['Lelm=',num2str(Lelms(5))],...
    ['Lelm=',num2str(Lelms(6))]
    };
lt = {'-', ':', ':', '--', '--','-'};

% length = 7;
% red = [1, 0, 0];
% pink = [0, 1, 1];
% colors_p = [linspace(red(1),pink(1),length)', linspace(red(2),pink(2),length)', linspace(red(3),pink(3),length)'];

% 	"#D95319"	
% Sample of RGB triplet [0.8500 0.3250 0.0980], which appears as dark orange
% 
% 	"#EDB120"	
% Sample of RGB triplet [0.9290 0.6940 0.1250], which appears as dark yellow
% 
% 	"#7E2F8E"	
% Sample of RGB triplet [0.4940 0.1840 0.5560], which appears as dark purple
% 
% 	"#77AC30"	
% Sample of RGB triplet [0.4660 0.6740 0.1880], which appears as medium green
% colors_p = [[0.4660 0.6740 0.1880]', ...
%             [0.9290 0.6940 0.1250]', ...
%             [0.3010 0.7450 0.9330]', ...
%             [0 0 1]', ...
%             [0.4940 0.1840 0.5560]', ...
%             [1 0 0]', ...
%             [0.8500 0.3250 0.0980]']'


figure_res = figure(1);
drawnow;
figure_res.Position = [200 100 750 400];
for iii = 1:6
    drawnow
    plot(x_axis_dmm, db(abs(LelmRes_Norm(iii,:))), 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii)));%, 'color', colors_p(iii,:))
    hold on
end
drawnow
plot(x_axis_dmm, db(abs(LelmRes_Norm(end,:))), '--','DisplayName',string(L(end)))%, 'color', colors_p(end,:))
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm]")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 5])
xlim([-15 15])
yticks(-60:6:6)
grid on

%%

TotResNorm = zeros(size(TotRes));
TotResNorm(:,:) = abs(TotRes(:,:));

[X,Y] = meshgrid(x_axis_dmm, 1:63);

% fig_surf = figure(2);
% fig_surf.Position = [200 100 1000 500];
% drawnow
% s = surf(X,Y,db(abs(TotResNorm)));
% s.EdgeColor = 'none'; 
% title("Lateral line through 30[mm]")
% xlabel("x [mm]")
% yticks(1:63)
% ylabel("Lelm")
% zlabel("Amplitude [dB]")
% xlim([-15 15])
% zlim([-60 0])


fig_2D = figure(3);
fig_2D.Position = [200 100 500 500];
drawnow
imagesc(x_axis_dmm,1:63,db(abs(TotResNorm)));
title("Lateral line through 30[mm]")
xlabel("x [mm]")
yticks([1 Lelms 62])
yticklabels({'DAS', '', num2str(Lelms(2)), num2str(Lelms(3)), num2str(Lelms(4)), num2str(Lelms(5)), num2str(Lelms(6))})
ylabel("Subarray size")
cb = colorbar();
% clim([-60 0]);
cb.Label.String = "Amplitude [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 4;
xlim([-10 10])

TotResNorm_all = zeros(size(TotRes));

for ii = 1:size(TotRes,1)
    TotResNorm_all(ii,:) = abs(TotRes(ii,:))./max(abs(TotRes(ii,:)));
end
fig_2D_norm = figure(4);
fig_2D_norm.Position = [800 100 500 500];
drawnow
imagesc(x_axis_dmm,1:63,db(abs(TotResNorm_all)));
title("Lateral line through 30[mm]")
xlabel("x [mm]")
yticks([1 Lelms 62])
yticklabels({'DAS', '', num2str(Lelms(2)), num2str(Lelms(3)), num2str(Lelms(4)), num2str(Lelms(5)), num2str(Lelms(6))})
ylabel("Subarray size")
cb = colorbar();
clim([-60 0]);
cb.Label.String = "Amplitude [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 4;
xlim([-10 10])

if save
    exportgraphics(figure_res, '..\Figures\Overleaf\Lelm_Resolution.pdf')
    exportgraphics(fig_2D, '..\Figures\Overleaf\Lelm_Resolution_2D.pdf')
    exportgraphics(fig_2D_norm, '..\Figures\Overleaf\Lelm_Resolution_2D_norm.pdf')
end


%% Figurer fra samme analyse med EV til tilsvarende subarraylengder
Lelms = 2:10:52;

fig_a_EV_Lelm = figure(5);
fig_a_EV_Lelm.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');

i = 1;
for Lelm = Lelms
    genname_az = "..\Figures\EV_Lelm\EV_az_"+num2str(Lelm)+"Lelm.fig";
    fig=openfig(genname_az);
    pause(0.5)

    xlim([25 34])
    ylim([1 10])
    colorbar off
    ax = gca;
    ax.Title.String = ['L = ', num2str(Lelm)];
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

    
    i = i+1;
    close(fig)

end
fig_a_EV_Lelm;
sgtitle("Eigenvalues of R through axial line at 0 degrees", 'FontSize', 18)

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Title.String = "[dB]";
% cb.Label.Rotation = 270;
% cb.Label.Position(1) = 4;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'Depth [mm]';
% tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
% tcl.YLabel.FontSize = 14;

fig_d_EV_Lelm = figure(6);
fig_d_EV_Lelm.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');

i = 1;
for Lelm = Lelms
    genname_az = "..\Figures\EV_Lelm\EV_dp_"+num2str(Lelm)+"Lelm.fig";
    fig=openfig(genname_az);
    pause(0.5)

    ylim([1 5])
    colorbar off
    ax = gca;
    ax.Title.String = ['L = ', num2str(Lelm)];
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

    i = i+1;
    close(fig)

end
fig_d_EV_Lelm;
sgtitle("Eigenvalues of R at lateral line through 30[mm]", 'FontSize', 18)

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Title.String = "[dB]";
% cb.Label.Rotation = 270;
% cb.Label.Position(1) = 4;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 16;


%% Ddist analyse filer
% Lfrac = 1/3 her, så omtrent 21 elementer per subarray i full lengde. K_in_lambda = 2
% Resten av verdiene er klassiske. 
fig_ddist_Amp = openfig(fullfile('..\Figures', '29-Mar-2023', 'Ddist', 'resolution_AmpCapon.fig'));
fig_ddist_Pow = openfig(fullfile('..\Figures', '29-Mar-2023', 'Ddist', 'resolution_PowCapon.fig'));
% tmpax = gca(fig_ddist_Pow);
% tmpax.YLabel.String = 'Amplitude';
fig_ddist_DAS = openfig(fullfile('..\Figures', '29-Mar-2023', 'Ddist', 'resolution_DAS.fig'));

dist = [0, 0.5, 1, 1.5, 2, 5];

fig_a_EV_dist = figure(100);
fig_a_EV_dist.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname_az      = "..\Figures\29-Mar-2023\Eigenvalues\Azim\EV_"  +num2str(d)  +   "dist_azim.fig";
    fig=openfig(genname_az);
    pause(0.5)

%     clim([-60 0])
    xlim([25 34])
    ylim([1 5])
%     cbb = colorbar()
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d)];
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

    
    i = i+1;
    close(fig)

end
fig_a_EV_dist;
sgtitle("Eigenvalues of R through axial line at 0 degrees", 'FontSize', 18)

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Title.String = "[dB]";
% cb.Label.String = "[dB]";
% cb.Label.Rotation = 270;
% cb.Label.Position(1) = 4;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'Depth [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 16;

fig_d_EV_dist = figure(101);
fig_d_EV_dist.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname_az      = "..\Figures\29-Mar-2023\Eigenvalues\Depth\EV_"  +num2str(d)  +   "dist_depth.fig";
    fig=openfig(genname_az);
    pause(0.5)

    ylim([1 5])
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d)];
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

    i = i+1;
    close(fig)

end
fig_d_EV_dist;
sgtitle("Eigenvalues of R at lateral line through 30[mm]", 'FontSize', 18)

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Title.String = "[dB]";

% cb.Label.String = "[dB]";
% cb.Label.Rotation = 270;
% cb.Label.Position(1) = 5;
cb.set('Limits', [-60 5])

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 16;

%% Condition Number ( CN cn cndl )
p = "C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\10-Mar-2023\Condition Number";
% CN_DL_0  = openfig(fullfile(p, "CNDL_1DL.fig"));
% % CN_DL_00001  = openfig(fullfile(p, "CNDL_0.0001DL.fig"));
% CN_DL_0001  = openfig(fullfile(p, "CNDL_0.001DL.fig"));
% CN_DL_001   = openfig(fullfile(p, "CNDL_0.01DL.fig"));
% CN_DL_01    = openfig(fullfile(p, "CNDL_0.1DL.fig"));
% CN_DL_1     = openfig(fullfile(p, "CNDL_1DL.fig"));

%% Plot med alle CN Max verdiene for en uendelig mengde DL verdier
pmax = "C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\FigFormat";
CN_DL_max  = openfig(fullfile(pmax, "maxCNDL.fig"));
% Må bruke lagrede data fra serverkjøring som går for øyeblikket
axObjs = CN_DL_max.Children;
dataObjs = axObjs.Children;
x = dataObjs.XData;
y = dataObjs.YData;
close(CN_DL_max)
% CN_DL_Max_log = figure();
% CN_DL_Max_log.Position = [100 100 700 350];
% plot(100*x, log10(y));
% xlabel("Diagonal load (%)")
% ylabel('log10(max(CN))')
% title("Maximum ConditionNumber for each diagonal load added to R")

path = "Data\06-Apr-2023";
filename = "CN_DL_Lelm_analysis.mat";

datamat = open(fullfile(path, filename));

CN_mid = datamat.CN_DL_Lelm_Mid;
CN_max = datamat.CN_DL_Lelm_Max;
DLs = datamat.DL;

tmpinds = find(CN_max > 1e30);
CN_max(tmpinds) = 0;
% CN_mid(tmpinds) = 0;

% totCNLelm = figure();
% subplot(121)
% imagesc(log10(CN_max))
% yticks(1:size(CN_max,1))
% yticklabels(split(num2str(DLs*100)))
% xlabel("Subarray size")
% ylabel("Diagonal load (%)")
% title("Maximum condition number in scene")
% c = colorbar();
% c.Label.String = "Condition number";
% 
% 
% % figure(2)
% subplot(122)
% imagesc(log10(CN_mid))
% yticks(1:size(CN_mid,1))
% yticklabels(split(num2str(DLs*100)))
% xlabel("Subarray size")
% ylabel("Diagonal load (%)")
% title("Condition number between scatterers in scene")
% c = colorbar();
% c.Label.String = "Condition number";

%% CNDLLelm linjeplot

Lelms = round([1/4, 1/3, 1/2, 2/3, 3/4]*64);
legends = {['L=', num2str(Lelms(1))], ['L=', num2str(Lelms(2))], ['L=', num2str(Lelms(3))], ['L=', num2str(Lelms(4))], ['L=', num2str(Lelms(5))]};
CNDLMAX = figure(1003);
CNDLMAX.Position = [50 100 750 400];

plot(DLs*100, log10(CN_max(:,Lelms)))

legend(legends, 'Location', 'northeast')
xlabel("Diagonal load (\%)")
ylabel('log10(CN)')
title("Maximum condition number of R in scene")
ylim([0 22])

breakyaxis([8 19]);

%%
CNDLMID = figure(1000);
CNDLMID.Position = [100 100 800 500];
subplot(3, 1, [1 2])
plot(DLs*100, log10(CN_mid(:,Lelms)))
legend(legends, 'Location', 'northeast')
% xlabel("Diagonal load (%)")
ylabel('log10(CN)')
xticks([0 5 10:10:100])
title("Condition Number of R between scatterers")

% CNDLMID_zoom = figure(1001);
% CNDLMID_zoom.Position = [50 100 700 350];
subplot(313)
plot(DLs*100, log10(CN_mid(:,Lelms)))
% legend(legends, 'Location', 'northeast')
xlabel("Diagonal load (\%)")
ylabel('log10(CN)')
% title("Condition Number of R between scatterers")
xlim([0 1])

% CNDLMID2 = figure(1002);
% CNDLMID2.Position = [800 400 700 350];
% plot(DLs*100, 1./CN_mid(:,Lelms))
% xlim([0 10])
% legend(legends, 'Location', 'northwest')
% xlabel("Diagonal load (%)")
% ylabel('1/CN')
% title("Condition Number of R between scatterers")


if save
    exportgraphics(CNDLMID, '..\Figures\Overleaf\CNDL_Lelms_Midpoint.pdf')
%     exportgraphics(CNDLMID_zoom, '..\Figures\Overleaf\CNDL_Lelms_Midpoint_zoomed.pdf')
%     exportgraphics(CNDLMID2, '..\Figures\Overleaf\CNDL_Lelms_Midpoint2.pdf')
    exportgraphics(CNDLMAX, '..\Figures\Overleaf\CNDL_Lelms_Maxval.pdf')
end


%% MLA koder ( mla )
p = "C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\MLA";
MLA_1  = openfig(fullfile(p, "result_1MLA.fig"));
MLA_1.Position = [100 100 450 400];
xlim([-5 5])

%%
MLA_3  = openfig(fullfile(p, "result_3MLA.fig"));
MLA_3.Position = [100 100 450 400];
xlim([-5 5])
%%
MLA_5  = openfig(fullfile(p, "result_5MLA.fig"));
MLA_5.Position = [100 100 450 400];
xlim([-5 5])
%%
p = "C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\MLA";
MLA_comp  = openfig(fullfile(p, "resolution_MLA.fig"));
pause(0.5)
ylim([-25 3])

%% SaveSection
if save
    path_ = "Overleaf";
    fig_path_ = "Bin";
    path = "..\Figures\Overleaf";

    exportgraphics(fig_a_EV_dist, fullfile(path, "EV_a_dist.pdf"))
    exportgraphics(fig_d_EV_dist, fullfile(path, "EV_d_dist.pdf"))
    exportgraphics(fig_a_EV_Lelm, fullfile(path, "EV_a_Lelm.pdf"))
    exportgraphics(fig_d_EV_Lelm, fullfile(path, "EV_d_Lelm.pdf"))

    exportgraphics(fig_ddist_DAS, fullfile(path, "res_ddist_DAS.pdf"))    
    exportgraphics(fig_ddist_Amp, fullfile(path, "res_ddist_Amp.pdf"))
    exportgraphics(fig_ddist_Pow, fullfile(path, "res_ddist_Pow.pdf"))

    p1 = "CNDL";
%     save_fig(CN_DL_0,     "CNDL_0.0DL",    type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")
% %     save_fig(CN_DL_00001, "CNDL_0.0001DL", type = ".pdf", path = fullpath("Overleaf", p1), fig_path="Bin")
%     save_fig(CN_DL_0001,  "CNDL_0.001DL",  type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")
%     save_fig(CN_DL_001,   "CNDL_0.01DL",   type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")
%     save_fig(CN_DL_01,    "CNDL_0.1DL",    type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")
%     save_fig(CN_DL_1,     "CNDL_1DL",      type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")

    p1 = "MLA";
    save_fig(MLA_1,    "result_1MLA",    type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")
    save_fig(MLA_3,    "result_3MLA",    type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")
    save_fig(MLA_5,    "result_5MLA",     type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")
    save_fig(MLA_comp, "resolution_MLA", type = ".pdf", path = fullfile("Overleaf", p1), fig_path="Bin")

end

%% Samling av figur scripts
% figs_speckle_levels
% % dataspeckle_FigureScript
