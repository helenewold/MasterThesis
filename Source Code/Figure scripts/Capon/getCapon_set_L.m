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
filename = 'Datasets/P4_FI_121444_45mm_focus.uff';

receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;
MLA = 1;
MLA_overlap = 1;

K_in_lambda = 0.5;
Lelm = 8;
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


%% getCapon USTB postprocess
% With set L, without normalization
Lelm_set=1;
USTB_normalization = 0;

script_post_getCapon
fig1 = b_data_mv_getCapon.plot(1, 'set L, no USTB normalization')
bdata_L_set_compare                =   uff.beamformed_data(b_data_mv_getCapon);
bdata_L_set_compare.data(:,1,1,1)  =   b_data_mv_getCapon.data      ./ max(b_data_mv_getCapon.data(:));
%% getCapon USTB postprocess
% With set L, with normalization
Lelm_set=1;
USTB_normalization = 1;

script_post_getCapon
fig2 = b_data_mv_getCapon.plot(2, 'set L, with USTB normalization');
bdata_L_set_compare.data(:,1,1,2)  =   b_data_mv_getCapon.data(:)    ./ max(b_data_mv_getCapon.data(:));


%% getCapon USTB postprocess
% Without set L, with normalization
Lelm_set=0;
USTB_normalization = 1;

script_post_getCapon

fig3 = b_data_mv_getCapon.plot(3, 'not set L, with USTB normalization');
bdata_L_set_compare.data(:,1,1,3)  =   b_data_mv_getCapon.data(:)    ./ max(b_data_mv_getCapon.data(:));
%% getCapon USTB postprocess
% Without set L, without normalization
Lelm_set=0;
USTB_normalization = 0;

script_post_getCapon

fig4 = b_data_mv_getCapon.plot(4, 'not set L, no USTB normalization');
bdata_L_set_compare.data(:,1,1,4)  =   b_data_mv_getCapon.data(:)    ./ max(b_data_mv_getCapon.data(:));

%% Movie loop comparison
% Movie loop of above plots.

% bdata_L_set_compare.plot(4,['1. Lset 2. Lset, norm 3. nothing 4. norm'])

% Notat her: ingen grunn til å normalisere slik når Lset = 1, fordi M vil
% ikke ha påvirkning på subarray length annet enn starten. Vil endre på men
% vil lagre resultater herfra for å kunne ta vare på dette. Må dokumentere
% denne endringen fordi Lelm_set kommer til å ekskludere ustb-normalisering
% uansett hva


%% 

savefig(fig2, "P4_FI_121444_45mm_focus_setL_normalized")
savefig(fig3, "P4_FI_121444_45mm_focus_dynL_normalized")
savefig(fig1, "P4_FI_121444_45mm_focus_setL")
savefig(fig4, "P4_FI_121444_45mm_focus_dynL")
