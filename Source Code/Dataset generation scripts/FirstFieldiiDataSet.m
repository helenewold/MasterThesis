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
filename = 'Field II dataset/MainTest3_degs.uff';
save = 1;
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

% %% USTB postprocess capon minimum variance
% script_post_capon_minvar

%% ustb getcapon postprocess
Lelm_set=0;
USTB_normalization = 1;
script_post_getCapon

% %% Plot Results
% % getCapon with apodization, plotted using imagesc, outside of USTB
% fig1 = figure(1);
% imagesc(azimuth_axis*1e3, depth_axis*1e3, db(abs(reshape(post_getCapon.data,post_getCapon.scan.N_depth_axis,post_getCapon.scan.N_azimuth_axis,post_getCapon.N_channels))))
% colormap('gray')
% title("getCapon with apodization");
% ylabel("z [mm]")
% xlabel("x [mm]")
% C = colorbar();
% C.Label.String = "Amplitude [dB]";
% 

% getCapon as USTB postprocess, using same apodization
fig2 = b_data_mv_getCapon.plot(2, 'getCapon postprocess');


% USTB capon_minimum_variance using same apodization
% fig3 = b_data_mv_MLA.plot(3, 'USTB Capon MV');

%% Movie loop of above plots.
% b_data_Capon_compare                =   uff.beamformed_data(b_data_mv_MLA);
% b_data_Capon_compare.data(:,1,1,1)  =   b_data_mv_MLA.data      ./ max(b_data_mv_MLA.data(:));
% b_data_Capon_compare.data(:,1,1,2)  =   post_getCapon.data(:)    ./ max(post_getCapon.data(:));
% % b_data_Capon_compare.data(:,1,1,3)  =   post_getCapon.data(:)    ./ max(post_getCapon.data(:));
% 
% fig4 = b_data_Capon_compare.plot(4,'1 = USTB Capon MV, 2 = getCapon postprocess');

%% Attempt of creating a diff-image. Calculation may be bad??
% b_data_diff = uff.beamformed_data(b_data_mv_MLA);
% data_diff = abs(b_data_mv_MLA.data(:) - post_getCapon.data(:));
% b_data_diff.data(:,1,1,1)  =   data_diff ./ max(data_diff(:));
% fig5 = b_data_diff.plot(5,'diff(USTB.getCapon, USTB.capon\_minvar)');

%% Figure saving block
if save
%     savefig(fig1, 'P4_FI_121444_45mm_focus_getCapon')
    savefig(fig2, 'FirstFieldiiDataSet')
%     savefig(fig3, 'P4_FI_121444_45mm_focus_USTB_capon_mv')
%     savefig(fig5, 'P4_FI_121444_45mm_focus_Capon_diffs')
end


