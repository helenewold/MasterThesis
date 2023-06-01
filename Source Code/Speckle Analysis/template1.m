% if ispc
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\downloads\Field_II_ver_3_30_windows.tar'));
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\Dokumenter\USTB'));
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\MasterThesis'));
% elseif isunix
%     addpath(genpath('/hom/dsb/field'));
%     addpath(genpath('/uio/hume/student-u55/helenewo/pc/Dokumenter/USTB'));
%     addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis'));
%     addpath(genpath('/uio/hume/studentg-u55/helenewo/MasterThesis/datasets'));
% end

%% Template to combine two datasets

% path = 'datasets';
% file_name = 'speckle_simulation_v1.uff';


filename = fullfile('Datasets', 'speckle_simulation_v1.uff');

% save_fig = 0;
% ow = 0;
% 
% plt = 0;


receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;
MLA = 1;
MLA_overlap = 1;

% Capon values (Kan endres underveis i scriptet)
K_in_lambda = 0.5;
L_frac = 1/3;
regCoeff = 1/100;

%% Read data
script_ReadData
%% Midprocess - Delay and sum (transmit dimension)
% b_data_MLA_delayed
script_mid_DAS_MLA
speckle_b_data = b_data_MLA_delayed;
%%
script_mid_DAS_bothDims
b_data_DAS_speckle = b_data_DAS;
fig1 = figure(1);
b_data_DAS_speckle.plot(fig1, "Simulated speckle");
%%
p = 30;%[20, 30, 40];
d = 2;
% dist = [0, 0.5, 1, 1.5, 2, 5];

degs = 360 * (d / (2*pi*p));

L_frac = 1/3;

MLA = 1;

phantom_positions(1,:) = [ p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
phantom_positions(2,:) = [-p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
        
h.save = 0;
h.amp = "ones";
    
channel_data = fieldii_generate_dataset("tmp", "tmp", phantom_positions,h);
%% Beamforming section - Kun getCapon beamforming
Lelm = channel_data.probe.N*L_frac;
%% Read data seksjon
% Setter altså variabler slik som i read data scriptet, men uten å lese fra
% fil.

channel_data.N_frames = 1;


script_mid_DAS_MLA
pkt_b_data = b_data_MLA_delayed;
script_mid_DAS_bothDims
b_data_DAS_pkt = b_data_DAS; 

b_data_DAS_pkt.plot(figure(2))
%% Combine data

b_data_MLA_delayed.data(:,:) = pkt_b_data.data(:,:) + 0.1*speckle_b_data.data(:,:);


%% ustb getcapon postprocess
Lelm_set=0;
calc_param = 0;
script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();
%%
b_data_mv_getCapon.plot(figure(3))
