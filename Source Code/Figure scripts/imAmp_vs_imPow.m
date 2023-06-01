%% Script explanation
% Skriptet viser resultatene som generelt viser hvordan få ut imPower fra
% getCapon postprocess. Tilsvarende endring må til i getCapon separat.
% Se script for kjøringen. 
% Datasettet er det som har vært brukt på starten, og parameterne er 
% uendret fra eksempelfilen gitt av Håvard som er brukt som utgangspunkt
% til USTB-prossesseringen.
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
% L_frac = 1/3;
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
Lelm_set=0;
USTB_normalization = 0;
script_post_getCapon

b_data_mv_getCapon_imAmplitude = mv_getCapon.go();

fig1 = b_data_mv_getCapon_imAmplitude.plot(1, 'imAmplitude');

%% Input imPower to object
b_data_mv_getCapon_imPower = uff.beamformed_data(b_data_mv_getCapon_imAmplitude);
b_data_mv_getCapon_imPower.data(:,:) = mv_getCapon.imPower(:,:);
fig2 = b_data_mv_getCapon_imPower.plot(2, 'imPower', 100);%, [0 100], [], [], [], [] , []);


%%
[imPower, Xs_p, Zs_p] = getScanConvertedImage(mv_getCapon.imPower, azimuth_axis, depth_axis, length(azimuth_axis), length(depth_axis));%, sizeX, sizeZ, interpolationMethod);
imPower = imPower./max(max(imPower));
[imAmplitude, Xs_a, Zs_a] = getScanConvertedImage(b_data_mv_getCapon_imAmplitude.data, azimuth_axis, depth_axis, length(azimuth_axis), length(depth_axis));%, sizeX, sizeZ, interpolationMethod);
imAmplitude = imAmplitude./max(max(imAmplitude));

figure(14)
subplot(121)
imagesc(Xs_p*1e3, Zs_p*1e3, db(abs(imPower)))
xlabel("Width [mm]")
ylabel("Depth [mm]")
title("imPower")
colormap(gray)
colorbar()

subplot(122)
imagesc(Xs_a*1e3, Zs_a*1e3, db(abs(imAmplitude)))
xlabel("Width [mm]")
ylabel("Depth [mm]")
title("imAmplitude")
colormap(gray)
colorbar()




%%

% savefig(fig1, "P4_FI_121444_45mm_focus_imAmp")
% savefig(fig2, "P4_FI_121444_45mm_focus_imPow")




