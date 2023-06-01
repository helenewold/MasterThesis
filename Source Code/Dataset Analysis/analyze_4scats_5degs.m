%% Parameter creation
path = 'Field II datasets';
file_name = 'FieldII_noSpeckle4scats_5degs.uff';

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

K_in_lambda = 2;
% Lelm = channel_data.probe.N/3;
L_frac = 1/3;
regCoeff = 1/100;

save_fig = 0;

%% Read data
script_ReadData

if ~exist('Lelm', 'var')
    Lelm = channel_data.probe.N*L_frac;
end
%% Steg to; midprocess - Delay and sum (NOT mla) (???)
script_mid_DAS_MLA

%% DAS begge dimensjoner - til bruk for Histogram Matching
% b_data_DAS
script_mid_DAS_bothDims
%% Steg tre: (optional) reshape vekk til 1 ramme
data_cube = USTB_reshape(scan_MLA, b_data_MLA_delayed, channel_data, n_frame = 1);

%% Danner getCapon postprocess structure
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
% 
% % getCapon som original prosess.
% K_samples = TimeAverageCalculation(scan_MLA, channel_data, K_in_lambda);
% script_getCapon;
% %
% getCapon_Scan = uff.beamformed_data(b_data_MLA_delayed);
% getCapon_Scan.data(:,1,1,1) = post_getCapon.data(:)./ max(post_getCapon.data(:));

%% USTB postprocess capon minimum variance
script_post_capon_minvar

%% ustb getcapon postprocess
Lelm_set=0;
calc_param = 0;
script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();

%% Histogram matching
[img_out, b_data_out] = tools.histogram_match(b_data_DAS,b_data_mv_getCapon);

%% Condition Number
if exist("calc_param", 'var') && calc_param == 1
    b_data_CN                   =   uff.beamformed_data(b_data_mv_getCapon);
    b_data_CN_DL             =   uff.beamformed_data(b_data_mv_getCapon);
    
    b_data_CN.data(:,:)         =   mv_getCapon.CN;
    b_data_CN_DL.data(:,:)   =   mv_getCapon.CN_DL;
end
%% Plotting section
% DAS, getCapon and USTB minimum variance results
fig1 = figure(1);
b_data_DAS.plot(fig1, "DAS output");

fig2 = figure(2);
b_data_mv_getCapon.plot(fig2, "getCapon output");

fig3 = figure(3);
b_data_mv_MLA.plot(fig3, "MV Capon output");

% Histogram matching

fig4 = figure(4);
b_data_out.plot(fig4, "Histogram matching (DAS - getCapon)", 60, 'none');
caxis([-60 0])

if exist("calc_param", 'var') && calc_param == 1

%% Condition number
fig5 = figure(5);
% fig5.WindowState = 'maximized';
% sgtitle('Condition numbers')
b_data_CN.plot(fig5, "Condition Number",[],'none');

fig6 = figure(6);
% fig6.WindowState = 'maximized';
sgtitle('Condition numbers - log10(data)')

subplot(121)
plot_CondNr(b_data_CN, azimuth_axis, depth_axis, title = "Before diagonal loading", dB='log10')

subplot(122)
plot_CondNr(b_data_CN_DL, azimuth_axis, depth_axis, title = "After diagonal loading", dB='log10')

fig7 = figure(7);
sgtitle('Condition numbers - percentage representation')
subplot(121)
plot_CondNr(b_data_CN, azimuth_axis, depth_axis, title = "Before diagonal loading", dB = '%')


subplot(122)
plot_CondNr(b_data_CN_DL, azimuth_axis, depth_axis, title = "After diagonal loading", dB = '%')
%% Percentage of eigenvalues above threshold
b_data_Reigval                  =   uff.beamformed_data(b_data_mv_getCapon);
b_data_Reigval.data(:,:)        =   mv_getCapon.Reigval;
b_data_Reigval_DL               =   uff.beamformed_data(b_data_mv_getCapon);
b_data_Reigval_DL.data(:,:)     =   mv_getCapon.Reigval_DL;

fig8 = figure(8);
b_data_Reigval.plot(fig8, '% of eigenvalues of R below 1e-10', [], 'none');

fig9 = figure(9);
b_data_Reigval_DL.plot(fig9, '% of eigenvalues of R below 1e-10, DL', [], 'none');

%% Tse values
b_data_Tse                  =   uff.beamformed_data(b_data_mv_getCapon);
b_data_Tse.data(:,:)        =   mv_getCapon.Tse;

fig10 = figure(10);
b_data_Tse.plot(fig10, 'T_{se} values', [], 'none');

%% White noise gain
b_data_WNG                  =   uff.beamformed_data(b_data_mv_getCapon);
b_data_WNG.data(:,:)        =   1./mv_getCapon.Tse;

fig11 = figure(11);
b_data_WNG.plot(fig11, 'White noise gain', [], 'none');



end