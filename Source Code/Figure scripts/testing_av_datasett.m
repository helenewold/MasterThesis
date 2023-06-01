%% NB
% Dette scriptet er bare for å kjøre randomly når jeg tester og ser om
% datasettet ble som ønmsket etter fieldii eller liknende. Ikke
% parameterendringer.

%% Setter parametre som skal brukes på tvers (Kan kanskje gjøres utenfor i eget script)
% filename = 'Datasets/P4_FI_121444_45mm_focus.uff';
filename = 'Source Code/Field II dataset/Final_Test3_degs.uff';

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
main_script

%% Plot Results
figure(1);
% sgtitle("getCapon")
% subplot(1,2,1)
imagesc(azimuth_axis*1e3, depth_axis*1e3, db(abs(reshape(post_getCapon.data,post_getCapon.scan.N_depth_axis,post_getCapon.scan.N_azimuth_axis,post_getCapon.N_channels))))
colormap('gray')
title("imagesc");
ylabel("z [mm]")
xlabel("x [mm]")
C = colorbar();
C.Label.String = "Amplitude [dB]";

% getCapon_Scan.plot(subplot(1,2,2), 'USTB plotting');


getCapon_Scan.plot(2, 'USTB plotting');

b_data_mv_MLA.plot(3, 'CAPON MLA');


b_data_Capon_compare                =   uff.beamformed_data(b_data_mv_MLA);
b_data_Capon_compare.data(:,1,1,1)  =   b_data_mv_MLA.data      ./ max(b_data_mv_MLA.data(:));
b_data_Capon_compare.data(:,1,1,2)  =   post_getCapon.data(:)    ./ max(post_getCapon.data(:));

b_data_Capon_compare.plot(4,['1 = USTB, 2 = getCapon'])

b_data_mv_getCapon.plot(5,['USTB.getCapon'])

