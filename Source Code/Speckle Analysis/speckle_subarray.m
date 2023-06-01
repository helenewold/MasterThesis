%% speckle_template
% Dette blir en template for hvordan regne ut allerede fullførte analyser
% innenfor subarray averaging mengde
close all
clear all
if isunix
    addpath(genpath('/hom/dsb/field'));
    addpath(genpath('/uio/hume/student-u55/helenewo/pc/Dokumenter/USTB'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/datasets'));
    path_ = ['./Figures/Server/',date];
    path_fig = ['./FiguresFigFormat/Server/',date];
    basepath = '/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis';
elseif ispc
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\downloads\Field_II_ver_3_30_windows.tar'));
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\Dokumenter\USTB'));
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\MasterThesis'));
% 
%     basepath = '\\hume.uio.no\student-u55\helenewo\MasterThesis\MasterThesis';

    basepath = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode';
end
% datapath = fullfile(basepath, 'Data', 'Fullrun', 'SpeckleLevel0.25', 'PartLelmAnalysis');
% EV_path = fullfile(basepath, 'Figures', 'Fullrun', 'SpeckleLevel0.25', 'PartLelmAnalysis');
% if ~exist(datapath); mkdir(datapath); end
% if ~exist(EV_path); mkdir(EV_path); end

% EV_path = fullfile(basepath, 'Figures', 'Fullrun', 'SpeckleLevel0.25', 'PartLelmAnalysis');
% if ~exist(EV_path); mkdir(EV_path); end

%% Load Data seksjon
% last inn speckle datasett
filename = fullfile('Datasets', 'speckle_simulation_v1.uff');
% save_fig = 0;
% ow = 1;
receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;

MLA = 3;
MLA_overlap = 1;

K_in_lambda = 2;

%% Read data
script_ReadData
channel_data_speckle = channel_data;
%%
p = 30;
depth_axis_full = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';
depth_inds = find(abs(depth_axis_full-p*1e-3)<1e-4);
depth_inds = depth_inds(round(length(depth_inds)/2));

depth_axis = depth_axis_full((depth_inds-40):(depth_inds+40));
%% DAS transmit dimension
% b_data_MLA_delayed
script_mid_DAS_MLA
speckle_b_data = b_data_MLA_delayed;
%% Generere 2pkt data

p = 30;%[20, 30, 40];
d = 2;%[0, 2, 5, 10, 20];
degs = 360 * (d / (2*pi*p));

% DL = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 1.5, 2, 2.5, 5:5:100]/100;
regCoeff = 1/100;

phantom_positions = zeros(2, 3);    
phantom_positions(1,:) = [ p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
phantom_positions(2,:) = [-p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
h.save = 0;
h.amp = "ones";

channel_data = fieldii_generate_dataset("tmp", "tmp", phantom_positions,h);

%% Beamforming section
Lelms = 42;%[2, 12, 22, 32, 42, 52];

channel_data.N_frames = 1;
script_mid_DAS_MLA

pkt_b_data = b_data_MLA_delayed;
b_data_MLA_delayed.data(:,:) = pkt_b_data.data(:,:) + 0.25*speckle_b_data.data(:,:);

cc = postprocess.coherent_compounding();
cc.input = b_data_MLA_delayed;
b_data_DAS = cc.go();
%%
% Lelms = [2, 12, 22, 32, 42, 52];%*1/2;%22:22;%2:63
lenL = length(Lelms);

LelmRes             = zeros(length(Lelms)+1, length(azimuth_axis)*MLA);

TotRes              = zeros(lenL+1, length(azimuth_axis)*MLA);
% TotRes_speckle      = zeros(lenL+1, length(azimuth_axis)*MLA);
% TotRes_3pkt         = zeros(lenL+1, length(azimuth_axis)*MLA);

LelmRes_Pow         = zeros(length(Lelms), length(azimuth_axis)*MLA);
TotRes_Pow          = zeros(lenL, length(azimuth_axis)*MLA);

CN_results          = zeros(8, lenL);

x_axis_dmm = p*azimuth_MLA;


%%
ii=2;
for i_Lelm = 1:length(Lelms)
    Lelm = Lelms(i_Lelm);
    disp(['Subarray length ', num2str(Lelm), ' of ', num2str(channel_data.N_channels), ' total'])
    
    %% ustb getcapon postprocess på fullt datasett
    Lelm_set=0;
    calc_param = 0;
    script_post_getCapon
    b_data_mv_getCapon = mv_getCapon.go();
    
    %% Plot eigenvalues
    K_samp = mv_getCapon.K_samples;
    N = mv_getCapon.scan.N_depth_axis;
    E = mv_getCapon.scan.N_azimuth_axis;
    az = deg2rad(0); % Rad
    az_inds = find(abs(azimuth_MLA-az)<0.01);
    az_inds = az_inds(round(length(az_inds)/2));
    
    depth_inds = find(abs(depth_axis-p*1e-3)<1e-4);
    depth_inds = depth_inds(round(length(depth_inds)/2));

    res_inds = depth_inds:length(depth_axis):length(b_data_mv_getCapon.data);


    max_length_az = max(cellfun('length',{mv_getCapon.StructTmp(:,az_inds).R_eigvals}));
    max_length_depth = max(cellfun('length',{mv_getCapon.StructTmp(depth_inds,:).R_eigvals}));
    Reig_axial = zeros(N, max_length_az);
    Reig_depth = zeros(E, max_length_depth);
    
    for i = 1:N
        len = cellfun('length',{mv_getCapon.StructTmp(i,az_inds).R_eigvals});
        Reig_axial(i, 1:len) = sort(mv_getCapon.StructTmp(i,az_inds).R_eigvals(:), 'descend');
        Reig_axial(i,(len+1):end) = nan;
    end
    for i = 1:E
        len = cellfun('length',{mv_getCapon.StructTmp(depth_inds,i).R_eigvals});
        Reig_depth(i, 1:len) = sort(mv_getCapon.StructTmp(depth_inds,i).R_eigvals(:), 'descend');
        Reig_depth(i,(len+1):end) = nan;
    end

    % Replace nan values
    Reig_axial(isnan(Reig_axial)) = -0.1;
    Reig_depth(isnan(Reig_depth)) = -0.1;

    fighandle1 = figure; 
    drawnow
    imagesc(depth_axis*1e3, 1:max_length_az, db(Reig_axial).');
    ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
    xlabel("Depth [mm]")
    title("Eigenvalues of R through " + num2str(az) + " degrees, Lelm="+num2str(Lelm))
    colbar = colorbar();
    colbar.Title.String = "Eigenvalue [dB]";

    fighandle2 = figure; 
    drawnow
    imagesc(x_axis_dmm, 1:max_length_depth, db(Reig_depth).')
    colbar = colorbar();
    colbar.Title.String = "Eigenvalue [dB]";
    title("Eigenvalues of R through 30 mm, Lelm="+num2str(Lelm))
    ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
    xlabel("Width [mm]")

%     savefig(fighandle1, fullfile(EV_path, ['EV_az_',num2str(Lelm),'Lelm']), 'compact')
%     savefig(fighandle2, fullfile(EV_path, ['EV_dp_',num2str(Lelm),'Lelm']), 'compact')
    %% Resolution through points
    tmp1 = depth_inds:length(depth_axis):length(b_data_mv_getCapon.data);
    tmp2 = (depth_inds-1):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp3 = (depth_inds-2):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp4 = (depth_inds+1):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp5 = (depth_inds+2):length(depth_axis):length(b_data_mv_getCapon.data);
    
    resolution_indexes = sort([tmp1 tmp2 tmp3 tmp4 tmp5]);
    for ind = 1:length(tmp1)
        tmpresind = sort([tmp1(ind) tmp2(ind) tmp3(ind) tmp4(ind) tmp5(ind)]);
        TotRes(Lelm, ind) = max((b_data_mv_getCapon.data(tmpresind)));
        TotRes_Pow(Lelm, ind) = max((mv_getCapon.imPower(tmpresind)));
    end

    if ismember(Lelm, Lelms)
        LelmRes(ii, :) = TotRes(Lelm, :);%./max(b_data_mv_getCapon.data(res_inds));
        LelmRes_Pow(ii, :) = TotRes_Pow(Lelm, :);%./max(b_data_mv_getCapon.data(res_inds));
        ii=ii+1;
    end
    %% Condition Number
    tmpaz = N*az_inds:1:N*(az_inds+1);
    tmpdepth = depth_inds:N:N*E;
    midpoint = intersect(tmpaz, tmpdepth);
    
    az_point = deg2rad(0.5*degs); % Rad
    az_i_R = find(abs(azimuth_MLA-az_point)<0.01);
    az_i_R = az_i_R(round(length(az_i_R)/2));
    rightpoint = intersect(N*az_i_R:1:N*(az_i_R+1), tmpdepth);
    az_i_L = find(abs(azimuth_MLA+az_point)<0.01);
    az_i_L = az_i_L(round(length(az_i_L)/2));
    leftpoint = intersect(N*az_i_L:1:N*(az_i_L+1), tmpdepth);

    CN = zeros(N,E);

    for i = 1:N
        for j = 1:E
            CN(i,j) = cond(mv_getCapon.StructTmp(i,j).R_DL);
        end
    end
%     CN_speckle = CN_speckle(:);
    b_data_CondNr_DL             =   uff.beamformed_data(b_data_mv_getCapon);
    b_data_CondNr_DL.data(:,:)   =   CN(:);

    CN_notnan = b_data_CondNr_DL.data(~isinf(b_data_CondNr_DL.data(~isnan(b_data_CondNr_DL.data))));
%     CN_speckle_notnan = CN_speckle(~isinf(CN_speckle(~isnan(CN_speckle))));
    % Maximum CN value storing
    CN_results(1, Lelm) = max(max(CN_notnan));
    % CN value between points storing
    CN_results(2, Lelm) = b_data_CondNr_DL.data(midpoint);
    % CN value On each point
    CN_results(3, Lelm) = b_data_CondNr_DL.data(leftpoint);
    CN_results(4, Lelm) = b_data_CondNr_DL.data(rightpoint);
    % CN value above middle
    CN_results(5, Lelm) = b_data_CondNr_DL.data(midpoint-8);
    CN_results(6, Lelm) = b_data_CondNr_DL.data(midpoint+8);


end
DAStmp = b_data_DAS.data(:);
tmp1 = depth_inds:length(depth_axis):length(DAStmp);
tmp2 = (depth_inds-1):length(depth_axis):length(DAStmp);
tmp3 = (depth_inds-2):length(depth_axis):length(DAStmp);
tmp4 = (depth_inds+1):length(depth_axis):length(DAStmp);
tmp5 = (depth_inds+2):length(depth_axis):length(DAStmp);
for ind = 1:length(tmp1)
    tmpresind = sort([tmp1(ind) tmp2(ind) tmp3(ind) tmp4(ind) tmp5(ind)]);
    TotRes(1, ind) = max(abs(DAStmp(tmpresind)));
end
LelmRes(1,:) = TotRes(1,:);
TotRes(1, :) = DAStmp(res_inds);

maxDas = max(abs(DAStmp(:)));

%%
% if ~exist(fullfile(basepath, 'Data', date)); mkdir(fullfile(basepath, 'Data', 'Fullrun')); end
% save(fullfile(datapath,'Lelm_analysis.mat'), 'TotRes', 'TotRes_Pow', 'LelmRes', 'LelmRes_Pow', 'maxDas', 'Lelms', 'CN_results', 'x_axis_dmm')
% save(fullfile(basepath, 'Data', 'Fullrun', 'Lelm_analysis.mat'), 'TotRes', 'TotRes_Pow', 'LelmRes', 'LelmRes_Pow', 'maxDas', 'Lelms', 'CN_results', 'x_axis_dmm')
% save(fullfile(basepath, 'Data', date, 'eigvals_subarray.mat'), 'eigmaxes', 'eigmins', 'eigmaxes_speckle', 'eigmins_speckle')