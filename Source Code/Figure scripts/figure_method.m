close all
clear all

fp = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Figures';
basepath = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode';
save = 1;
%% prosessering av data for å lage figurene at all
p = 30;
d = 2; degs = 360 * (d / (2*pi*p));

Lfrac = 1/3;

MLA = 3;

phantom_positions(1,:) = [ p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
phantom_positions(2,:) = [-p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
    
h.save = 0;
h.amp = "ones";

channel_data = fieldii_generate_dataset("tmp", "tmp", phantom_positions,h);

receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;
MLA_overlap = 1;

regCoeff = 1/100;
K_in_lambda = 2;

% Definerer scan
azimuth_axis=zeros(channel_data.N_waves,1);
depth_axis = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';

for n=1:channel_data.N_waves
    azimuth_axis(n)=channel_data.sequence(n).source.azimuth;
end  

Lelm = channel_data.probe.N*Lfrac;

channel_data.N_frames = 1;

script_mid_DAS_MLA
b_data_pkt_DASRec = b_data_MLA_delayed;

cc = postprocess.coherent_compounding();
cc.input = b_data_pkt_DASRec;
b_data_DAS_pkt = cc.go();
%%
filename = fullfile('Datasets', 'speckle_simulation_v1.uff');
script_ReadData

script_mid_DAS_MLA
b_data_speckle_DASRec = b_data_MLA_delayed;

cc = postprocess.coherent_compounding();
cc.input = b_data_speckle_DASRec;
b_data_DAS_speckle = cc.go();
%%
b_data_MLA_delayed.data(:,:) = b_data_pkt_DASRec.data(:,:) + 0.25*b_data_speckle_DASRec.data(:,:);
b_data_full_DASRec = b_data_MLA_delayed;

cc = postprocess.coherent_compounding();
cc.input = b_data_full_DASRec;
b_data_DAS_full = cc.go();
%%
MLA = 1;
filename = fullfile('Datasets', 'phantom_cyst.uff');
script_ReadData

script_mid_DAS_bothDims
b_data_cyst_DAS = b_data_DAS;


%%
MLA = 1;
filename = fullfile('Datasets', 'P4_FI_121444_45mm_focus.uff');
script_ReadData

script_mid_DAS_bothDims
b_data_phantom_DAS = b_data_DAS;

%% Plotting
fig_pkt = figure(1);
fig_speckle = figure(2);
fig_full = figure(3);
fig_cyst = figure(4);
fig_phantom = figure(5);
b_data_DAS_pkt.plot(fig_pkt, 'Simulated dataset with d=2mm between scatterers')
b_data_DAS_speckle.plot(fig_speckle, 'Simulated speckle dataset')
b_data_DAS_full.plot(fig_full, 'Full simulated dataset')
b_data_cyst_DAS.plot(fig_cyst, 'Simulated cyst dataset')
b_data_phantom_DAS.plot(fig_phantom, 'Phantom dataset')

%% 
if save
    exportgraphics(fig_pkt,     fullfile(fp, 'Overleaf', 'Method', 'dataset_point.pdf'))
    exportgraphics(fig_speckle, fullfile(fp, 'Overleaf', 'Method', 'dataset_speckle.pdf'))
    exportgraphics(fig_full,    fullfile(fp, 'Overleaf', 'Method', 'dataset_full.pdf'))
    exportgraphics(fig_cyst,    fullfile(fp, 'Overleaf', 'Method', 'dataset_cyst.pdf'))
    exportgraphics(fig_phantom, fullfile(fp, 'Overleaf', 'Method', 'dataset_phantom.pdf'))
end

