%% Script explanation
% Skriptet viser resultatene som oppsummerer arbeidet gjort med å få
% getCapon over i USTB postprocess, og å legge på apodisering både i
% postprosessen og i getCapon separat. Se script for kjøringen. Datasettet
% er det som har vært brukt på starten, og parameterne er uendret fra
% eksempelfilen gitt av Håvard som er brukt som utgangspunkt til
% USTB-prossesseringen.
%
% Oppdatert: 26.10.2022

%% Parameter creation
path = 'Datasets';
file_name = 'P4_FI_121444_45mm_focus';
filename = append(fullfile(path, file_name),'.uff');

save = 0;


receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;
MLA = 1;
MLA_overlap = 1;

K_in_lambda = 0.5;
% Lelm = channel_data.probe.N/3;
L_frac = 1/3;
regCoeff = 1/100;


%% Read data
script_ReadData

%% Steg to; midprocess - Delay and sum (NOT mla) (???)
script_mid_DAS_MLA

%% Steg tre: (optional) reshape vekk til 1 ramme
data_cube = USTB_reshape(scan_MLA, b_data_MLA_delayed, channel_data, n_frame = 1);

%% Danner getCapon postprocess structure
if ~exist('Lelm', 'var')
    Lelm = channel_data.probe.N*L_frac;
end

post_getCapon.dimension = dimension.receive;

post_getCapon.transmit_apodization= mid.transmit_apodization;
post_getCapon.receive_apodization = mid.receive_apodization;
post_getCapon.scan = scan_MLA;

post_getCapon.channel_data = channel_data;

post_getCapon.K_in_lambda = K_in_lambda;
post_getCapon.L_elements = Lelm;
post_getCapon.regCoef = regCoeff;

post_getCapon.input = b_data_MLA_delayed;

%% Apodization matrix
[rx_apodization, tx_apodization] = create_apod_matrix(post_getCapon);

%% getCapon som original prosess.
K_samples = TimeAverageCalculation(scan_MLA, channel_data, K_in_lambda);
script_getCapon;

getCapon_Scan = uff.beamformed_data(b_data_MLA_delayed);
getCapon_Scan.data(:,1,1,1) = post_getCapon.data(:)./ max(post_getCapon.data(:));

%% USTB postprocess capon minimum variance
script_post_capon_minvar
%% ustb getcapon postprocess
Lelm_set=0;
USTB_normalization = 1;
calc_param = 0;

script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();

%%
subarrsize_ustb = mv_MLA.subarrsize;
subarrsize_getCapon = mv_getCapon.subarrsize;
%%
fig123 = figure(123);
% subplot(131)
% imagesc(subarrsize_ustb(10:502,:))
% subplot(132)
% imagesc(subarrsize_getCapon(10:502,:))
% subplot(133)
imagesc(azimuth_axis, depth_axis*1e3, subarrsize_ustb(10:502,:)-subarrsize_getCapon(10:502,:))
ylabel("z[mm]")
xlabel("x [rad]")
colormap('gray')
colorbar

savefig("C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\getCaponLelmDiff.fig")
exportgraphics(fig123, "C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\Overleaf\getCaponLelmDiff.pdf")
%% Plot Results
% getCapon with apodization, plotted using imagesc, outside of USTB
fig1 = figure(1);
imagesc(azimuth_axis*1e3, depth_axis*1e3, db(abs(reshape(post_getCapon.data,post_getCapon.scan.N_depth_axis,post_getCapon.scan.N_azimuth_axis,post_getCapon.N_channels))))
colormap('gray')
title("getCapon with apodization");
ylabel("z [mm]")
xlabel("x [mm]")
C = colorbar();
C.Label.String = "Amplitude [dB]";


% getCapon as USTB postprocess, using same apodization
fig2 = figure(2);
b_data_mv_getCapon.plot(fig2, 'getCapon postprocess');


% USTB capon_minimum_variance using same apodization
fig3 = figure(3);
b_data_mv_MLA.plot(fig3, 'USTB Capon MV');


%% Movie loop of above plots.
b_data_Capon_compare                =   uff.beamformed_data(b_data_mv_MLA);
b_data_Capon_compare.data(:,1,1,1)  =   b_data_mv_MLA.data      ./ max(b_data_mv_MLA.data(:));
b_data_Capon_compare.data(:,1,1,2)  =   post_getCapon.data(:)    ./ max(post_getCapon.data(:));
% b_data_Capon_compare.data(:,1,1,3)  =   post_getCapon.data(:)    ./ max(post_getCapon.data(:));

fig4 = figure(4);
b_data_Capon_compare.plot(fig4,'1 = USTB Capon MV, 2 = getCapon postprocess');

%% Attempt of creating a diff-image. Calculation may be bad??
b_data_diff = uff.beamformed_data(b_data_mv_MLA);
diff1 = abs(b_data_mv_MLA.data./max(b_data_mv_MLA.data(:)));
diff2 = abs(post_getCapon.data./max(post_getCapon.data(:)));

% data_diff = abs(b_data_mv_MLA.data(:) - post_getCapon.data(:));
b_data_diff.data(:,1,1,1)  = abs(diff1(:)-diff2(:));

fig5 = figure(5);
b_data_diff.plot(fig5,'Difference between getCapon and capon\_minvar');
c = colorbar();


clear diff1 diff2

%% Figure saving block
if exist('save','var') && save
    savefig(fig1, "getCapon_"       + file_name,     path = "Capon_Section\"+file_name)
    savefig(fig2, "getCaponPost_"   + file_name,     path = "Capon_Section\"+file_name)
    savefig(fig3, "MVcapon_"        + file_name,     path = "Capon_Section\"+file_name)
    savefig(fig5, "diff_"           + file_name,     path = "Capon_Section\"+file_name)
end

%%
figtst = openfig("C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\getCaponLelmDiff.fig")
mapcolors = [0 0 0; 1 1 1];
cmap = [repmat(mapcolors(1,:),[1 1]); ...
    mapcolors(2,:)];
set(gca,'clim',[0 1])
colormap(cmap);
c = colorbar;
c.XTickLabel = {'', '0', '', '1', ''}
c.XTick = [0:1/4:1]

exportgraphics(figtst, "C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\Overleaf\getCaponLelmDiff.pdf")

