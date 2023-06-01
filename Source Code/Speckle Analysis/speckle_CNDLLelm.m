%% speckle_template_CNDL
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
    addpath(genpath('\\hume.uio.no\studeaznt-u55\helenewo\pc\downloads\Field_II_ver_3_30_windows.tar'));
    addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\Dokumenter\USTB'));
    addpath(genpath('\\hume.uio.no\student-u55\helenewo\MasterThesis'));
    basepath = '\\hume.uio.no\student-u55\helenewo\MasterThesis\MasterThesis';
end
% if ~exist('speckle_value')
speckle_value = 0.25;
% end
data_path = fullfile(basepath, 'Data', 'Fullrun', ['SpeckleLevel', num2str(speckle_value)], 'CNDLLelm');
if ~exist(data_path); mkdir(data_path); end

%% Load Data seksjon
filename = fullfile('datasets', 'speckle_simulation_v1.uff');
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

% Read data
script_ReadData
channel_data_speckle = channel_data;

p = 30;
depth_axis_full = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';
depth_axis_30mmInds = find(abs(depth_axis_full-p*1e-3)<1e-4);
depth_axis_30mmInds = depth_axis_30mmInds(round(length(depth_axis_30mmInds)/2));

depth_axis = depth_axis_full((depth_axis_30mmInds-50):(depth_axis_30mmInds+50));

% DAS transmit dimension
% b_data_MLA_delayed
script_mid_DAS_MLA
speckle_b_data = b_data_MLA_delayed;

%% 2pt
p = 30;%[20, 30, 40];
d = 2;%[0, 2, 5, 10, 20];
degs = 0.5 * 360 * (d / (2*pi*p));

DL = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 1.5, 2, 2.5, 5:5:100]/100;
% DL = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5]%, 1, 1.5, 2, 2.5, 5:5:100]/100;
% DL = 1/100;
CN_max = zeros(length(DL), 1);
phantom_positions = zeros(2, 3);

phantom_positions(1,:) = [ p*1e-3*sind( degs ), 0, p*1e-3*cosd( degs )];
phantom_positions(2,:) = [-p*1e-3*sind( degs ), 0, p*1e-3*cosd( degs )];
            
h.save = 0;
h.amp = "ones";
        
channel_data = fieldii_generate_dataset("tmp", "tmp", phantom_positions,h);

script_mid_DAS_MLA
pkt_b_data = b_data_MLA_delayed;

b_data_MLA_delayed.data(:,:) = pkt_b_data.data(:,:) + speckle_value*speckle_b_data.data(:,:);


%% Beamforming section - Kun getCapon beamforming

channel_data.N_frames = 1;

cc = postprocess.coherent_compounding();
cc.input = b_data_MLA_delayed;
b_data_DAS = cc.go();

x_axis_dmm = p*(azimuth_MLA);

CN_results = zeros(8, length(DL), length(1:63));
% CN_Speckle_results = zeros(8, length(DL), length(1:63));
eigmaxes            = zeros(length(DL), 63, length(azimuth_axis)*MLA*length(depth_axis));
eigmins             = zeros(length(DL), 63, length(azimuth_axis)*MLA*length(depth_axis));

TotRes = zeros(length(DL),63, length(azimuth_axis)*MLA);
TotRes_Pow = zeros(length(DL),63, length(azimuth_axis)*MLA);
LDL = length(DL)
assignin('base', 'LDL', LDL)

%%
depth_axis_30mmInds_new = find(abs(depth_axis-p*1e-3)<1e-4);
depth_axis_30mmInds_new = depth_axis_30mmInds_new(round(length(depth_axis_30mmInds_new)/2));

N = length(depth_axis);
E = length(azimuth_axis)*MLA;
az = deg2rad(0); % Rad

az_inds = find(abs(azimuth_MLA-az)<0.01);
az_inds = az_inds(round(length(az_inds)/2));

% depth_axis_30mmInds = find(abs(depth_axis-p*1e-3)<1e-4);
% depth_axis_30mmInds = depth_axis_30mmInds(round(length(depth_axis_30mmInds)/2));
tmpaz = N*az_inds:1:N*(az_inds+1);
tmpdepth = depth_axis_30mmInds_new:N:N*E;
midpoint = intersect(tmpaz, tmpdepth);


az_point = deg2rad(0.5*degs); % Rad
az_i_R = find(abs(azimuth_MLA-az_point)<0.01);
az_i_R = az_i_R(round(length(az_i_R)/2));
rightpoint = intersect(N*az_i_R:1:N*(az_i_R+1), tmpdepth);
az_i_L = find(abs(azimuth_MLA+az_point)<0.01);
az_i_L = az_i_L(round(length(az_i_L)/2));
leftpoint = intersect(N*az_i_L:1:N*(az_i_L+1), tmpdepth);


tmp1 = depth_axis_30mmInds_new:N:N*E;
tmp2 = (depth_axis_30mmInds_new-1):N:N*E;
tmp3 = (depth_axis_30mmInds_new-2):N:N*E;
tmp4 = (depth_axis_30mmInds_new+1):N:N*E;
tmp5 = (depth_axis_30mmInds_new+2):N:N*E;
   
%%
p = gcp('nocreate'); % If no pool, do not create new one.

if isempty(p)
    myCluster=parcluster('local');
    myCluster.NumWorkers=45;
    parpool(myCluster,myCluster.NumWorkers)

end
%%

parfor DL_i = 1:LDL
    disp(['Work started on worker no ', num2str(DL_i)])
    regCoeff = DL(DL_i);
    
    for Lelm = 2:63
        Lelm_set=0;

        calc_param = 0;
        mv_getCapon = postprocess.capon_MV_getCapon();
        mv_getCapon.dimension = dimension.receive;
        
        mv_getCapon.transmit_apodization= mid.transmit_apodization;
        mv_getCapon.receive_apodization = mid.receive_apodization;
        mv_getCapon.scan = scan_MLA;
        
        mv_getCapon.channel_data = channel_data;
        
        mv_getCapon.K_in_lambda = K_in_lambda;
        mv_getCapon.L_elements = Lelm;
        

        mv_getCapon.L_elements_set = 0;
        mv_getCapon.USTB_normalization = 1;
        mv_getCapon.calc_param = 0;        
        
        mv_getCapon.regCoef = regCoeff;
        
        mv_getCapon.input = b_data_MLA_delayed;
        
        [rx_apodization, tx_apodization] = create_apod_matrix(mv_getCapon);
        
        b_data_mv_getCapon = mv_getCapon.go();

        %%
        CN = zeros(N,E);
        for i = 1:N
            for j = 1:E
                CN(i,j) = cond(mv_getCapon.StructTmp(i,j).R_DL);
            end
        end
        b_data_CondNr_DL             =   uff.beamformed_data(b_data_mv_getCapon);
        b_data_CondNr_DL.data(:,:)   =   CN(:);

        CN_notnan = b_data_CondNr_DL.data(~isinf(b_data_CondNr_DL.data(~isnan(b_data_CondNr_DL.data))));

        CN_results(:,DL_i, Lelm) = [max(max(CN_notnan)), b_data_CondNr_DL.data(midpoint), b_data_CondNr_DL.data(leftpoint), b_data_CondNr_DL.data(rightpoint), b_data_CondNr_DL.data(midpoint-8), b_data_CondNr_DL.data(midpoint+8), b_data_CondNr_DL.data(N*E*1/4), b_data_CondNr_DL.data(N*E*3/4)].';

        %% Resolution through points

        for ind = 1:E
            tmpresind = sort([tmp1(ind) tmp2(ind) tmp3(ind) tmp4(ind) tmp5(ind)]);
            TotRes(DL_i, Lelm, ind) = max((b_data_mv_getCapon.data(tmpresind)));
            TotRes_Pow(DL_i, Lelm, ind) = max((mv_getCapon.imPower(tmpresind)));
        end
        disp(['Lelm no. ', num2str(Lelm), ' of 64 done calculating at DL no. ', num2str(DL_i), ' of ', num2str(length(DL))])
    end
    disp(['DL no. ', num2str(DL_i), ' of ', num2str(length(DL)), ' done calculating'])
end

DasRes = zeros(length(azimuth_axis)*MLA, 1);
for ind_tmp = 1:E
    tmpresind = sort([tmp1(ind_tmp) tmp2(ind_tmp) tmp3(ind_tmp) tmp4(ind_tmp) tmp5(ind_tmp)]);
    DasRes(ind_tmp) = max(b_data_DAS.data(tmpresind));
end
maxDas = max(abs(b_data_DAS.data(:)));
save(fullfile(data_path, 'CN_DL_Lelm_speckle.mat'), 'x_axis_dmm', 'maxDas', 'CN_results', 'DL', 'TotRes', 'TotRes_Pow', "DasRes")











