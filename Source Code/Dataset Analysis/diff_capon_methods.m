close all
clear


%% Parameter creation
path = 'Datasets';
file_name = 'speckle_simulation_v1'

filename = append(fullfile(path, file_name),'.uff');

save = 1;



receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;
MLA = 1;
MLA_overlap = 1;

% Capon values (Kan endres underveis i scriptet)
K_in_lambda = 2;
L_frac = 1/3;
regCoeff = 1/100;

%% Read data
script_ReadData

%% Midprocess - Delay and sum (transmit dimension)
% b_data_MLA_delayed
script_mid_DAS_MLA

%% (optional) reshape vekk til 1 frame
data_cube = USTB_reshape(scan_MLA, b_data_MLA_delayed, channel_data, n_frame = 1);

Lelm = channel_data.probe.N*L_frac;

%% USTB Capon
% b_data_mv_MLA
script_post_capon_minvar

%% getCapon postprocess
script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();


%% Diff image
b_data_diff = uff.beamformed_data(b_data_mv_MLA);

% data_diff = abs(b_data_mv_MLA.data(:) - post_getCapon.data(:));
b_data_diff.data(:,1,1,1)  = abs(abs(b_data_mv_MLA.data./max(b_data_mv_MLA.data(:))) ...
                        -abs(b_data_mv_getCapon.data./max(b_data_mv_getCapon.data(:))));
fig = figure(1);
b_data_diff.plot(fig, 'Difference between getCapon and USTB Capon')

if save
    exportgraphics(fig, 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\Overleaf\CaponDiffs.pdf')
    savefig(fig, 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures\Overleaf\CaponDiffs.fig')
end
