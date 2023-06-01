%% Script explanation
% En initiell template for histogram matching.

clear
close all
%% Parameter creation
path = 'Field II datasets';
file_name = 'FieldII_noSpeckle10scats_5degs.uff';

% path = 'Datasets';
% file_name = 'P4_FI_121444_45mm_focus.uff';
filename = fullfile(path, file_name);

receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;
MLA = 1;
MLA_overlap = 1;

K_in_lambda = 5;
% Lelm = channel_data.probe.N/3;
L_frac = 1/3;
regCoeff = 1/100;

save_fig = 0;

%% Read data
script_ReadData

%% Steg to; midprocess - Delay and sum
% b_data_MLA_delayed
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

%% USTB postprocess capon minimum variance
% b_data_mv_MLA
script_post_capon_minvar

%% ustb getcapon postprocess

script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();

%% DAS helt ut
% b_data_DAS
script_mid_DAS_bothDims

%% Plotte DAS resultat
fig1 = figure(1);
b_data_DAS.plot(fig1, "DAS output")

fig2 = figure(2);
b_data_mv_getCapon.plot(fig2, "getCapon output");

fig3 = figure(3);
b_data_mv_MLA.plot(fig3, "MV Capon output");

%% Histogram matching

[img_out, b_data_out] = tools.histogram_match(b_data_DAS,b_data_mv_getCapon);

fig4 = figure(4);
b_data_out.plot(fig4, "USTB plotting - Histogram matching", 60, 'none');
caxis([-60 0])

%% Figure saving block
if exist('save_fig','var') && save_fig
    savefig(fig1, "HM_DAS_"     + file_name,     path = "Histogram_Matching")
    savefig(fig4, "HM_USTB_"    + file_name,     path = "Histogram_Matching")
end


