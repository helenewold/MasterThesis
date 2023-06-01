clear;
close all;

path = "..\Field II Datasets";

degs_separated = -5;

pos = 5:5:50;

n_scats = length(pos);

phantom_positions = zeros(n_scats*2, 3);

for i = 1:n_scats
    phantom_positions(i,:) = [pos(i)/1000*sind( degs_separated ), 0, pos(i)/1000*cosd( degs_separated )];
    phantom_positions(n_scats+i,:) = [-pos(i)/1000*sind( degs_separated ), 0, pos(i)/1000*cosd( degs_separated )];
end

h.degs_separated = degs_separated;
h.save = 1;
h.amp = "ones";

filename = "FieldII_noSpeckle" + num2str(n_scats) + "scats_" + num2str(abs(degs_separated)) + "degs";
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


%% Midprocess - Delay and sum
script_mid_DAS_MLA

%% (optional) reshape vekk til 1 frame
data_cube = USTB_reshape(scan_MLA, b_data_MLA_delayed, channel_data, n_frame = 1);

%% getCapon postprocess structure
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

% Apodization matrix
[rx_apodization, tx_apodization] = create_apod_matrix(post_getCapon);



%% Postprocess - getCapon
script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();

%% Plot result
fig1 = figure(1);
b_data_mv_getCapon.plot(fig1, 'getCapon postprocess');
