%% generate_2pkt_ddist.m
% Script to generate points at a distance d from each other.
clear;
close all;

path = "..\Field II datasets";
filename = "d_dist_3x2points";

pos = [20, 30, 40];
d = 10;
n_scats = length(pos);

phantom_positions = zeros(n_scats*2, 3);

for i = 1:n_scats
    phantom_positions(i,:) = [0.5*d*1e-3, 0, pos(i)*1e-3];
    phantom_positions(n_scats+i,:) = [-0.5*d*1e-3, 0, pos(i)*1e-3];
end

h.save = 1;
h.amp = "ones";

channel_data = fieldii_generate_dataset(path, filename, phantom_positions,h);


%% Beamforming section - Kun getCapon beamforming
% Parameter creation
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

if ~exist('Lelm', 'var')
    Lelm = channel_data.probe.N*L_frac;
end

%% Read data seksjon
% Setter altså variabler slik som i read data scriptet, men uten å lese fra
% fil.

channel_data.N_frames = 1;

% Definerer scan
azimuth_axis=zeros(channel_data.N_waves,1);
depth_axis = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';

for n=1:channel_data.N_waves
    azimuth_axis(n)=channel_data.sequence(n).source.azimuth;
end
scan=uff.linear_scan('azimuth_axis',azimuth_axis,'depth_axis',depth_axis);


script_mid_DAS_MLA

%% DAS begge dimensjoner - til bruk for Histogram Matching
% b_data_DAS
script_mid_DAS_bothDims


%% USTB postprocess capon minimum variance
script_post_capon_minvar

%% ustb getcapon postprocess
Lelm_set=0;
calc_param = 0;
script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();

%% Histogram matching
[img_out, b_data_out] = tools.histogram_match(b_data_DAS,b_data_mv_getCapon);


%% Plotting
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


