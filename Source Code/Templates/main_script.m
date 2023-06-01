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
%%
getCapon_Scan = uff.beamformed_data(b_data_MLA_delayed);
getCapon_Scan.data(:,1,1,1) = post_getCapon.data(:)./ max(post_getCapon.data(:));

%% USTB postprocess capon minimum variance
script_post_capon_minvar

%% ustb getcapon postprocess
Lelm_set=0;
script_post_getCapon



