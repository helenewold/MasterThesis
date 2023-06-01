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
% elseif ispc
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\downloads\Field_II_ver_3_30_windows.tar'));
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\Dokumenter\USTB'));
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\MasterThesis'));
%     basepath = '\\hume.uio.no\student-u55\helenewo\MasterThesis\MasterThesis';
end
% if ~exist('speckle_value')
basepath = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode';
speckle_value = 0.25;
data_path = fullfile(basepath, 'Data', 'Contrast');
if ~exist(data_path); mkdir(data_path); end

%% Load Data seksjon
% filepath = fullfile('uio','hume','student-u55','helenewo','MasterThesis', 'datasets');
filename = fullfile('Datasets', 'phantom_cyst.uff');

% save_fig = 0;
% ow = 1;
receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;

MLA = 1;

MLA_overlap = 1;

K_in_lambda = 2;

% Read data
script_ReadData

% depth_axis = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';
%%
depth_axis_full = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';
depth_axis_30mmInds = find(abs(depth_axis_full-30*1e-3)<1e-4);
depth_axis_30mmInds = depth_axis_30mmInds(round(length(depth_axis_30mmInds)/2));

depth_axis = depth_axis_full%((depth_axis_30mmInds-90):(depth_axis_30mmInds+90));

% DAS transmit dimension
% b_data_MLA_delayed
script_mid_DAS_MLA
%% Finne koordinater

N = length(depth_axis);
E = length(azimuth_axis)*MLA;


%         ROI_coord = ((dybden*1e-3).^2)

xc_ROI = -0;
zc_ROI = 30;
r_ROI = 3.5;
r_ROI_small = 2.5;
r_background_inner = 6;
r_background_outer = 9;
rad_circ = 5;

z_position = zeros(N, E);
for azim = 1:E
    z_position(:, azim) = depth_axis(:);
end    
z_position = scan_MLA.z;%z_position(:);
x_position = zeros(N, E);
for depth = 1:N
    x_position(depth, :) = depth_axis(depth)*azimuth_MLA(:);
end
x_position = scan_MLA.x;%(x_position(:));

points = ((x_position - xc_ROI*1e-3).^2 + (z_position - zc_ROI*1e-3).^2);
idx_ROI_orig = points < (r_ROI*1e-3)^2;
idx_ROI_small = points < (r_ROI_small*1e-3)^2;
idx_ROI_circ = points < (rad_circ*1e-3)^2;


idx_background_outer =  (((x_position-xc_ROI*10^-3).^2) + (z_position-zc_ROI*10^-3).^2 < (r_background_outer*10^-3)^2); %ROI speckle
idx_background_inner =  (((x_position-xc_ROI*10^-3).^2) + (z_position-zc_ROI*10^-3).^2 < (r_background_inner*10^-3)^2); %ROI speckle
idx_background_outer(idx_background_inner) = 0;
idx_background = idx_background_outer;
    


%% Beamforming section - Kun getCapon beamforming

DL = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 1.5, 2, 2.5, 5:5:100]/100;
Lelms = 2:63;
%%
% DL = 1/100;
% Lelms = [10, 20, 40, 60];
% Lelms = 64/2

channel_data.N_frames = 1;

cc = postprocess.coherent_compounding();
cc.input = b_data_MLA_delayed;
b_data_DAS = cc.go();
%%
angs = 0:0.1:2*pi;
xa_circle       = rad_circ*sin(angs);
za_circle       = rad_circ*cos(angs) + zc_ROI;
xa_small       = r_ROI_small*sin(angs);
za_small       = r_ROI_small*cos(angs) + zc_ROI;
xa_inner_ROI    = r_ROI*sin(angs);
za_inner_ROI    = r_ROI*cos(angs) + zc_ROI;

xa_inner_B      = r_background_inner*sin(angs);
za_inner_B      = r_background_inner*cos(angs) + zc_ROI;
xa_outer_B      = r_background_outer*sin(angs);
za_outer_B      = r_background_outer*cos(angs) + zc_ROI;
fig1 = figure(1);
b_data_DAS.plot(fig1, 'DAS beamformed data')
plot(xa_circle, za_circle, 'r:', 'LineWidth',2)
plot(xa_small, za_small, 'g:', 'LineWidth',2)
plot(xa_inner_ROI, za_inner_ROI, 'y:', 'LineWidth',2)
plot(xa_inner_B, za_inner_B, 'b:', 'LineWidth',2)
plot(xa_outer_B, za_outer_B, 'b:', 'LineWidth',2)
exportgraphics(fig1, fullfile(basepath, 'Figures', 'Overleaf', 'Contrast', 'result_marked.pdf'))
savefig(fig1, fullfile(basepath, 'Figures', 'Overleaf', 'Contrast', 'result_marked.fig'))

%%
gCNR_orig     = zeros(length(DL),63);
CNR_orig      = zeros(length(DL),63);
CR_orig       = zeros(length(DL),63);

gCNR_circ     = zeros(length(DL),63);
CNR_circ      = zeros(length(DL),63);
CR_circ       = zeros(length(DL),63);

gCNR_small     = zeros(length(DL),63);
CNR_small      = zeros(length(DL),63);
CR_small       = zeros(length(DL),63);


LDL = length(DL);
assignin('base', 'LDL', LDL);

  
%%
p = gcp('nocreate'); % If no pool, do not create new one.

if isempty(p)
    myCluster=parcluster('local');
    myCluster.NumWorkers=45;
    parpool(myCluster,myCluster.NumWorkers)

end
%%
gcp
parfor DL_i = 1:LDL
    disp(['Work started on worker no ', num2str(DL_i)])
    regCoeff = DL(DL_i);
    
    for Lelm = Lelms
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
%         fig = figure();
%         b_data_mv_getCapon.plot(fig, [num2str(Lelm), ' elements'])
%         exportgraphics(fig, fullfile(basepath, 'Figures', 'Overleaf', 'Contrast', ['result_', num2str(Lelm), 'Lelms.pdf']))
%         savefig(fig, fullfile(basepath, 'Figures', 'Overleaf', 'Contrast', ['result_', num2str(Lelm), 'Lelms.fig']))

        %% Kommenter ut for å kunne fikse lokalt
%         if Lelm == round(30) && regCoeff == 1/100
%             angs = 0:0.1:2*pi;
%             xa_circle       = rad_circ*sin(angs);
%             za_circle       = rad_circ*cos(angs) + zc_ROI;
%             xa_inner_ROI    = r_ROI*sin(angs);
%             za_inner_ROI    = r_ROI*cos(angs) + zc_ROI;
%             xa_inner_B      = r_background_inner*sin(angs);
%             za_inner_B      = r_background_inner*cos(angs) + zc_ROI;
%             xa_outer_B      = r_background_outer*sin(angs);
%             za_outer_B      = r_background_outer*cos(angs) + zc_ROI;
% 
%             tmp_background = uff.beamformed_data(b_data_mv_getCapon);
%             tmp_background.data = (idx_background + idx_ROI).*b_data_mv_getCapon.data;
% %             tmp_background.plot(figure(3))
%             fig1 = figure();
%             b_data_mv_getCapon.plot(fig1, "Cyst phantom, ROI and background sones marked")
%             plot(xa_circle, za_circle, 'r:', 'LineWidth',1.25)
%             plot(xa_inner_ROI, za_inner_ROI, 'y:', 'LineWidth',1.25)
%             plot(xa_inner_B, za_inner_B, 'b:', 'LineWidth',1.25)
%             plot(xa_outer_B, za_outer_B, 'b:', 'LineWidth',1.25)
%             exportgraphics(fig1, fullfile(basepath, 'Figures', 'Overleaf', 'Contrast', 'result_marked.pdf'))
%             savefig(fig1, fullfile(basepath, 'Figures', 'Overleaf', 'Contrast', 'result_marked.fig'))
% 
%         end

        %% Regne verdier
        data_background = b_data_mv_getCapon.data(idx_background);

        data_ROI_orig = b_data_mv_getCapon.data(idx_ROI_orig);
        gCNR_orig(DL_i, Lelm)     = gCNR(abs(data_ROI_orig), abs(data_background));
        CR_orig(DL_i, Lelm)       = CR(abs(data_ROI_orig), abs(data_background));       % output dB?
        CNR_orig(DL_i, Lelm)      = CNR(abs(data_ROI_orig), abs(data_background));

        data_ROI_circ = b_data_mv_getCapon.data(idx_ROI_circ);
        gCNR_circ(DL_i, Lelm)     = gCNR(abs(data_ROI_circ), abs(data_background));
        CR_circ(DL_i, Lelm)       = CR(abs(data_ROI_circ), abs(data_background));       % output dB?
        CNR_circ(DL_i, Lelm)      = CNR(abs(data_ROI_circ), abs(data_background));

        data_ROI_small = b_data_mv_getCapon.data(idx_ROI_small);
        gCNR_small(DL_i, Lelm)     = gCNR(abs(data_ROI_small), abs(data_background));
        CR_small(DL_i, Lelm)       = CR(abs(data_ROI_small), abs(data_background));       % output dB?
        CNR_small(DL_i, Lelm)      = CNR(abs(data_ROI_small), abs(data_background));


        %% Resolution through points

        disp(['Lelm no. ', num2str(Lelm), ' of 64 done calculating at DL no. ', num2str(DL_i), ' of ', num2str(length(DL))])
    end
    disp(['DL no. ', num2str(DL_i), ' of ', num2str(length(DL)), ' done calculating'])
end
%%
data_background_DAS =   b_data_DAS.data(idx_background);
data_ROI_DAS_orig   =   b_data_DAS.data(idx_ROI_orig);
data_ROI_DAS_circ   =   b_data_DAS.data(idx_ROI_circ);
data_ROI_DAS_small  =   b_data_DAS.data(idx_ROI_small);

gCNR_DAS_orig       =   gCNR(abs(data_ROI_DAS_orig), abs(data_background_DAS));
CR_DAS_orig         =   CR(abs(data_ROI_DAS_orig), abs(data_background_DAS));
CNR_DAS_orig        =   CNR(abs(data_ROI_DAS_orig), abs(data_background_DAS));

gCNR_DAS_circ       =   gCNR(abs(data_ROI_DAS_circ), abs(data_background_DAS));
CR_DAS_circ         =   CR(abs(data_ROI_DAS_circ), abs(data_background_DAS));
CNR_DAS_circ        =   CNR(abs(data_ROI_DAS_circ), abs(data_background_DAS));

gCNR_DAS_small      =   gCNR(abs(data_ROI_DAS_small), abs(data_background_DAS));
CR_DAS_small        =   CR(abs(data_ROI_DAS_small), abs(data_background_DAS));
CNR_DAS_small       =   CNR(abs(data_ROI_DAS_small), abs(data_background_DAS));



%%
save(fullfile(data_path,'contrast_metrics.mat'), ...
    'CNR_orig',      'CR_orig',      'gCNR_orig', ...
    'CNR_circ',      'CR_circ',      'gCNR_circ', ...
    'CNR_small',     'CR_small',     'gCNR_small', ...
    'CNR_DAS_orig',  'CR_DAS_orig',  'gCNR_DAS_orig', ...
    'CNR_DAS_circ' , 'CR_DAS_circ',  'gCNR_DAS_circ', ...
    'CNR_DAS_small', 'CR_DAS_small', 'gCNR_DAS_small', ...
    'DL', 'Lelms')











