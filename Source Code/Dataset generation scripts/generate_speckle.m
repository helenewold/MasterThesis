
clear;
close all;

path = "..\Field II datasets\Speckle";
filename = "scat_3kscats";


h.save = 1;



N = 3000;

x_size = 50/1000;   %  Width of phantom [mm]
y_size = 0/1000;   %  Transverse width of phantom [mm]
z_size = 58/1000;   %  Height of phantom [mm]
z_start = 15/1000;  %  Start of phantom surface [mm];

%  Create the general scatterers

x = (rand(N,1)-0.5)*x_size;
y = zeros(N,1);
z = rand(N,1)*z_size + z_start;
h.amp = "rand";
h.D = 32;
h.fs=100e6;
phantom_positions(:,:,:) = [x, y, z];

%%
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
%%
depth_testing = 177:512:65536;
tst = uff.beamformed_data(b_data_mv_getCapon);
tst.data(depth_testing) = 1;


fig5 = figure(5);
tst.plot(fig5, "getCapon output testing");
% hold on
% plot(20e-3*ones(128,1), depth_testing)
% hold off

