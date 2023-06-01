%% Script explanation
% Skriptet viser første resultat med en initiell bruk av Condition Number.
%
% Oppdatert: 10.11.22
%           11.11.22:   Endret plottemetode
%           21.11.22:   Oppdatert plotting og seksjoner
close all
clear all

%% Parameter creation
% path = 'Field II datasets';
% file_name = 'FieldII_noSpeckle4scats_5degs';

path = 'Datasets';
file_name = 'P4_FI_121444_45mm_focus';

filename = append(fullfile(path, file_name),'.uff');

save_fig = 0;


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


%% Read data
script_ReadData

%% Midprocess - Delay and sum
script_mid_DAS_MLA

%% (optional) reshape vekk til 1 frame
data_cube = USTB_reshape(scan_MLA, b_data_MLA_delayed, channel_data, n_frame = 1);

%% getCapon postprocess structure
if ~exist('Lelm', 'var')
    Lelm = channel_data.probe.N*L_frac;
end

% post_getCapon.dimension = dimension.receive;
% 
% post_getCapon.transmit_apodization= mid.transmit_apodization;
% post_getCapon.receive_apodization = mid.receive_apodization;
% post_getCapon.scan = scan_MLA;
% 
% post_getCapon.channel_data = channel_data;
% 
% post_getCapon.K_in_lambda = K_in_lambda;
% post_getCapon.L_elements = Lelm;
% post_getCapon.regCoef = regCoeff;
% 
% post_getCapon.input = b_data_MLA_delayed;
% 
% % Apodization matrix
% [rx_apodization, tx_apodization] = create_apod_matrix(post_getCapon);


%% Postprocess - getCapon
calc_param = 1;
script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();

%% Plot result
fig1 = figure(1);
b_data_mv_getCapon.plot(fig1, 'getCapon postprocess');

%% Condition number plotting
b_data_CondNr                   =   uff.beamformed_data(b_data_mv_getCapon);
b_data_CondNr_DL             =   uff.beamformed_data(b_data_mv_getCapon);

b_data_CondNr.data(:,:)         =   mv_getCapon.CN;
b_data_CondNr_DL.data(:,:)   =   mv_getCapon.CN_DL;

%%
fig2 = figure(2);
fig2.WindowState = 'maximized';
sgtitle('Condition numbers')

subplot(121)
% b_data_CondNr_DL.plot(subplot(121), 'tittel', [], 'none')
plot_CondNr(b_data_CondNr, azimuth_axis, depth_axis, title = "Before diagonal loading")%, dB = 'none')


subplot(122)
plot_CondNr(b_data_CondNr_DL, azimuth_axis, depth_axis, title = "After diagonal loading")%, dB = 'none')


%%
fig3 = figure(3);
fig3.WindowState = 'maximized';
sgtitle('Condition numbers - log10(data)')

subplot(121)
plot_CondNr(b_data_CondNr, azimuth_axis, depth_axis, title = "Before diagonal loading", dB='log10')

subplot(122)
plot_CondNr(b_data_CondNr_DL, azimuth_axis, depth_axis, title = "After diagonal loading", dB='log10')



%% Figure saving block
if exist('save_fig','var') && save_fig
    savefig(fig2, "CN_"     + file_name,     path = "Condition_Number")
end
