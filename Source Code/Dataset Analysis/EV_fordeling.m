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
else
    basepath = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode';
end


% EV_path = fullfile(basepath, 'Figures', date, 'LelmEV');
% if ~exist(EV_path); mkdir(EV_path); end

%% Load Data seksjon
% last inn speckle datasett
filename = fullfile('datasets', 'speckle_simulation_v1.uff');
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

%% Read data
script_ReadData
channel_data_speckle = channel_data;
%% DAS transmit dimension
% b_data_MLA_delayed
script_mid_DAS_MLA
speckle_b_data             =   uff.beamformed_data(b_data_MLA_delayed);

%% Generere 2pkt data

p = 30;%[20, 30, 40];
d = 2;%[0, 2, 5, 10, 20];
degs = 360 * (d / (2*pi*p));

% DL = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 1.5, 2, 2.5, 5:5:100]/100;
regCoeff = 1/100;

phantom_positions = zeros(3, 3);    
phantom_positions(1,:) = [ p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
phantom_positions(2,:) = [ 0, 0, p*1e-3*cosd( degs*0.5 )];
phantom_positions(3,:) = [-p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
h.save = 0;
h.amp = "ones";

channel_data = fieldii_generate_dataset("tmp", "tmp", phantom_positions,h);

%% Beamforming section
Lelms = [2, 12, 22, 32, 42, 52];

channel_data.N_frames = 1;
script_mid_DAS_MLA
% script_mid_DAS_bothDims
% pkt_b_data = b_data_MLA_delayed;

pkt_b_data             =   uff.beamformed_data(b_data_MLA_delayed);

b_data_MLA_delayed.data(:,:) = pkt_b_data.data(:,:) + 0.1*speckle_b_data.data(:,:);

cc = postprocess.coherent_compounding();
cc.input = b_data_MLA_delayed;
b_data_DAS = cc.go();
b_data_DAS.plot()

cc_3pkt = postprocess.coherent_compounding();
cc_3pkt.input = pkt_b_data;
b_data_DAS_3pkt = cc_3pkt.go();

cc_speckle = postprocess.coherent_compounding();
cc_speckle.input = speckle_b_data;
b_data_DAS_speckle = cc_speckle.go();

x_axis_dmm = p*azimuth_MLA;


%%
partLelms = round([1/10, 1/4, 1/3, 1/2, 2/3, 3/4]*64);%*1/2;%22:22;%2:63
lenL = length(partLelms);

LelmRes             = zeros(length(Lelms)+1, length(azimuth_axis)*MLA);

TotRes              = zeros(lenL+1, length(azimuth_axis)*MLA);
TotRes_speckle      = zeros(lenL+1, length(azimuth_axis)*MLA);
TotRes_3pkt         = zeros(lenL+1, length(azimuth_axis)*MLA);

LelmRes_Pow         = zeros(length(Lelms), length(azimuth_axis)*MLA);
TotRes_Pow          = zeros(lenL, length(azimuth_axis)*MLA);

CN_results          = zeros(8, lenL);
CN_3pkt_results     = zeros(8, lenL);
CN_Speckle_results  = zeros(8, lenL);

eigmins             = zeros(lenL, length(azimuth_axis)*MLA*length(depth_axis));
eigmaxes            = zeros(lenL, length(azimuth_axis)*MLA*length(depth_axis));
eigmins_3pkt        = zeros(lenL, length(azimuth_axis)*MLA*length(depth_axis));
eigmaxes_3pkt       = zeros(lenL, length(azimuth_axis)*MLA*length(depth_axis));
eigmins_speckle     = zeros(lenL, length(azimuth_axis)*MLA*length(depth_axis));
eigmaxes_speckle    = zeros(lenL, length(azimuth_axis)*MLA*length(depth_axis));

histcounts_tot_foc      = zeros(lenL, 100);
histcounts_tot_nonfoc   = zeros(lenL, 100);


%%
ii=2;
for i_Lelm = 1:length(partLelms)
    Lelm = partLelms(i_Lelm);
    disp(['Subarray length ', num2str(Lelm), ' of ', num2str(channel_data.N_channels), ' total'])
    
%     %% ustb getcapon postprocess på fullt datasett
%     Lelm_set=0;
%     calc_param = 0;
%     script_post_getCapon
%     b_data_mv_getCapon = mv_getCapon.go();
%     b_data_mv_getCapon.plot()
        %% USTB getCapon på kun speckle
    mv_getCapon_speckle = postprocess.capon_MV_getCapon();
    mv_getCapon_speckle.dimension = dimension.receive;
    
    mv_getCapon_speckle.transmit_apodization= mid.transmit_apodization;
    mv_getCapon_speckle.receive_apodization = mid.receive_apodization;
    mv_getCapon_speckle.scan = scan_MLA;
    
    mv_getCapon_speckle.channel_data = channel_data;
    
    mv_getCapon_speckle.K_in_lambda = K_in_lambda;
    mv_getCapon_speckle.L_elements = Lelm;
    
    
    if exist("Lelm_set", 'var')
        mv_getCapon_speckle.L_elements_set = Lelm_set;
    else
        mv_getCapon_speckle.L_elements_set = 0;
    end
    
    
    if exist("USTB_normalization", 'var')
        mv_getCapon_speckle.USTB_normalization = USTB_normalization;
    else
        mv_getCapon_speckle.USTB_normalization = 1;
    end
    
    if exist("calc_param", 'var')
        mv_getCapon_speckle.calc_param = calc_param;
    % else
    %     mv_getCapon_speckle.calc_param = 0;
    end
    
    % mv_getCapon_speckle.fig_handle = fighandle;
   
    
    mv_getCapon_speckle.regCoef = regCoeff;
    
    mv_getCapon_speckle.input = speckle_b_data;
    
    [rx_apodization, tx_apodization] = create_apod_matrix(mv_getCapon_speckle);

    b_data_mv_getCapon_Speckle = mv_getCapon_speckle.go();
        %% USTB getCapon på kun datapkt
%     mv_getCapon_3pkt = postprocess.capon_MV_getCapon();
%     mv_getCapon_3pkt.dimension = dimension.receive;
%     
%     mv_getCapon_3pkt.transmit_apodization= mid.transmit_apodization;
%     mv_getCapon_3pkt.receive_apodization = mid.receive_apodization;
%     mv_getCapon_3pkt.scan = scan_MLA;
%     
%     mv_getCapon_3pkt.channel_data = channel_data;
%     
%     mv_getCapon_3pkt.K_in_lambda = K_in_lambda;
%     mv_getCapon_3pkt.L_elements = Lelm;
%     
%     
%     if exist("Lelm_set", 'var')
%         mv_getCapon_3pkt.L_elements_set = Lelm_set;
%     else
%         mv_getCapon_3pkt.L_elements_set = 0;
%     end
%     
%     
%     if exist("USTB_normalization", 'var')
%         mv_getCapon_3pkt.USTB_normalization = USTB_normalization;
%     else
%         mv_getCapon_3pkt.USTB_normalization = 1;
%     end
%     
%     if exist("calc_param", 'var')
%         mv_getCapon_3pkt.calc_param = calc_param;
%     end
%     
%    
%     
%     mv_getCapon_3pkt.regCoef = regCoeff;
%     
%     mv_getCapon_3pkt.input = pkt_b_data;
%     
%     [rx_apodization, tx_apodization] = create_apod_matrix(mv_getCapon_3pkt);
% 
%     b_data_mv_getCapon_3pkt = mv_getCapon_3pkt.go();
    
%     %% Plot eigenvalues
%     K_samp = mv_getCapon.K_samples;
    N = mv_getCapon_speckle.scan.N_depth_axis;
    E = mv_getCapon_speckle.scan.N_azimuth_axis;
    az = deg2rad(0); % Rad
    az_inds = find(abs(azimuth_MLA-az)<0.01);
    az_inds = az_inds(round(length(az_inds)/2));
    
    depth_inds = find(abs(depth_axis-p*1e-3)<1e-4);
    depth_inds = depth_inds(round(length(depth_inds)/2));
% 
%     res_inds = depth_inds:length(depth_axis):length(b_data_mv_getCapon.data);
% 
% 
%     max_length_az = max(cellfun('length',{mv_getCapon.StructTmp(:,az_inds).R_eigvals}));
%     max_length_depth = max(cellfun('length',{mv_getCapon.StructTmp(depth_inds,:).R_eigvals}));
%     Reig_axial = zeros(N, max_length_az);
%     Reig_depth = zeros(E, max_length_depth);
%     
%     for i = 1:N
%         len = cellfun('length',{mv_getCapon.StructTmp(i,az_inds).R_eigvals});
%         Reig_axial(i, 1:len) = sort(mv_getCapon.StructTmp(i,az_inds).R_eigvals(:), 'descend');
%         Reig_axial(i,(len+1):end) = nan;
%     end
%     for i = 1:E
%         len = cellfun('length',{mv_getCapon.StructTmp(depth_inds,i).R_eigvals});
%         Reig_depth(i, 1:len) = sort(mv_getCapon.StructTmp(depth_inds,i).R_eigvals(:), 'descend');
%         Reig_depth(i,(len+1):end) = nan;
%     end
% 
%     % Replace nan values
%     Reig_axial(isnan(Reig_axial)) = -0.1;
%     Reig_depth(isnan(Reig_depth)) = -0.1;
% 
%     fighandle1 = figure; 
%     drawnow
%     imagesc(depth_axis*1e3, 1:max_length_az, db(Reig_axial).');
%     ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
%     xlabel("Depth [mm]")
%     title("Eigenvalues of R through " + num2str(az) + " degrees, Lelm="+num2str(Lelm))
%     colbar = colorbar();
%     colbar.Title.String = "Eigenvalue [dB]";
% 
%     fighandle2 = figure; 
%     drawnow
%     imagesc(x_axis_dmm, 1:max_length_depth, db(Reig_depth).')
%     colbar = colorbar();
%     colbar.Title.String = "Eigenvalue [dB]";
%     title("Eigenvalues of R through 30 mm, Lelm="+num2str(Lelm))
%     ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
%     xlabel("Width [mm]")

%     savefig(fighandle1, fullfile(EV_path, ['EV_az_',num2str(Lelm),'Lelm']), 'compact')
%     savefig(fighandle2, fullfile(EV_path, ['EV_dp_',num2str(Lelm),'Lelm']), 'compact')
    %%
%     tmp1 = depth_inds:length(depth_axis):length(b_data_mv_getCapon.data);
%     tmp2 = (depth_inds-1):length(depth_axis):length(b_data_mv_getCapon.data);
%     tmp3 = (depth_inds-2):length(depth_axis):length(b_data_mv_getCapon.data);
%     tmp4 = (depth_inds+1):length(depth_axis):length(b_data_mv_getCapon.data);
%     tmp5 = (depth_inds+2):length(depth_axis):length(b_data_mv_getCapon.data);
%     
%     resolution_indexes = sort([tmp1 tmp2 tmp3 tmp4 tmp5]);
%     for ind = 1:length(tmp1)
%         tmpresind = sort([tmp1(ind) tmp2(ind) tmp3(ind) tmp4(ind) tmp5(ind)]);
%         TotRes(i_Lelm+1, ind) = max((b_data_mv_getCapon.data(tmpresind)));
%         TotRes_3pkt(i_Lelm+1, ind) = max((b_data_mv_getCapon_3pkt.data(tmpresind)));
%         TotRes_speckle(i_Lelm+1, ind) = max((b_data_mv_getCapon_Speckle.data(tmpresind)));
%         TotRes_Pow(i_Lelm+1, ind) = max((mv_getCapon.imPower(tmpresind)));
%     end
% 
%     if ismember(Lelm, Lelms)
%         LelmRes(ii, :) = TotRes(i_Lelm+1, :);%./max(b_data_mv_getCapon.data(res_inds));
%         LelmRes_Pow(ii, :) = TotRes_Pow(i_Lelm+1, :);%./max(b_data_mv_getCapon.data(res_inds));
%         ii=ii+1;
%     end
    %% Condition Number
%     tmpaz = N*az_inds:1:N*(az_inds+1);
%     tmpdepth = depth_inds:N:N*E;
%     midpoint = intersect(tmpaz, tmpdepth);
%     
%     az_point = deg2rad(0.5*degs); % Rad
%     az_i_R = find(abs(azimuth_MLA-az_point)<0.01);
%     az_i_R = az_i_R(round(length(az_i_R)/2));
%     rightpoint = intersect(N*az_i_R:1:N*(az_i_R+1), tmpdepth);
%     az_i_L = find(abs(azimuth_MLA+az_point)<0.01);
%     az_i_L = az_i_L(round(length(az_i_L)/2));
%     leftpoint = intersect(N*az_i_L:1:N*(az_i_L+1), tmpdepth);
% 
%     CN = zeros(N,E);
%     CN_speckle = zeros(N,E);
%     CN_pkt = zeros(N,E);
%     eigmax = zeros(N,E);
%     eigmin = zeros(N,E);
%     eigmax_speckle = zeros(N,E);
%     eigmin_speckle = zeros(N,E);
%     eigmax_3pkt = zeros(N,E);
%     eigmin_3pkt = zeros(N,E);
%     for i = 1:N
%         for j = 1:E
%             CN(i,j) = cond(mv_getCapon.StructTmp(i,j).R_DL);
%             CN_speckle(i,j) = cond(mv_getCapon_speckle.StructTmp(i,j).R_DL);
%             CN_pkt(i,j) = cond(mv_getCapon_3pkt.StructTmp(i,j).R_DL);
% 
%             if ~isempty(mv_getCapon.StructTmp(i,j).R_eigvals_DL)
%                 eigmax(i,j) = max(mv_getCapon.StructTmp(i,j).R_eigvals_DL);
%                 eigmin(i,j) = min(mv_getCapon.StructTmp(i,j).R_eigvals_DL);
%             end
%             if ~isempty(mv_getCapon_speckle.StructTmp(i,j).R_eigvals_DL)
%                 eigmax_speckle(i,j) = max(mv_getCapon_speckle.StructTmp(i,j).R_eigvals_DL);
%                 eigmin_speckle(i,j) = min(mv_getCapon_speckle.StructTmp(i,j).R_eigvals_DL);
%             end
%             if ~isempty(mv_getCapon_3pkt.StructTmp(i,j).R_eigvals_DL)
%                 eigmax_3pkt(i,j) = max(mv_getCapon_3pkt.StructTmp(i,j).R_eigvals_DL);
%                 eigmin_3pkt(i,j) = min(mv_getCapon_3pkt.StructTmp(i,j).R_eigvals_DL);
%             end
%         end
%     end
%     CN_speckle = CN_speckle(:);
%     CN_pkt = CN_speckle(:);
%     b_data_CondNr_DL             =   uff.beamformed_data(b_data_mv_getCapon);
%     b_data_CondNr_DL.data(:,:)   =   CN(:);
% 
%     CN_notnan = b_data_CondNr_DL.data(~isinf(b_data_CondNr_DL.data(~isnan(b_data_CondNr_DL.data))));
%     CN_speckle_notnan = CN_speckle(~isinf(CN_speckle(~isnan(CN_speckle))));
%     % Maximum CN value storing
%     CN_results(1, i_Lelm) = max(max(CN_notnan));
%     % CN value between points storing
%     CN_results(2, i_Lelm) = b_data_CondNr_DL.data(midpoint);
%     % CN value On each point
%     CN_results(3,  i_Lelm) = b_data_CondNr_DL.data(leftpoint);
%     CN_results(4,  i_Lelm) = b_data_CondNr_DL.data(rightpoint);
%     % CN value above middle
%     CN_results(5,  i_Lelm) = b_data_CondNr_DL.data(midpoint-8);
%     CN_results(6,  i_Lelm) = b_data_CondNr_DL.data(midpoint+8);
%     % Random above_below
%     CN_results(7,  i_Lelm) = CN_notnan(N/4*E/4);
%     CN_results(8,  i_Lelm) = CN_notnan(N*3/4*E*3/4);
% 
%     CN_Speckle_results(1, i_Lelm) = max(max(CN_speckle_notnan));
%     CN_Speckle_results(2, i_Lelm) = CN_speckle(midpoint);
%     CN_Speckle_results(3, i_Lelm) = CN_speckle(leftpoint);
%     CN_Speckle_results(4, i_Lelm) = CN_speckle(rightpoint);
%     CN_Speckle_results(5, i_Lelm) = CN_speckle(midpoint-5);
%     CN_Speckle_results(6, i_Lelm) = CN_speckle(midpoint+5);
%     CN_Speckle_results(7, i_Lelm) = CN_speckle_notnan(N/4*E/4);
%     CN_Speckle_results(8, i_Lelm) = CN_speckle_notnan(N*3/4*E*3/4);
% 
%     CN_3pkt_results(1, i_Lelm) = max(max(CN_speckle_notnan));
%     CN_3pkt_results(2, i_Lelm) = CN_speckle(midpoint);
%     CN_3pkt_results(3, i_Lelm) = CN_speckle(leftpoint);
%     CN_3pkt_results(4, i_Lelm) = CN_speckle(rightpoint);
%     CN_3pkt_results(5, i_Lelm) = CN_speckle(midpoint-5);
%     CN_3pkt_results(6, i_Lelm) = CN_speckle(midpoint+5);
%     CN_3pkt_results(7, i_Lelm) = CN_speckle_notnan(N/4*E/4);
%     CN_3pkt_results(8, i_Lelm) = CN_speckle_notnan(N*3/4*E*3/4);
% 
%     eigmaxes(i_Lelm,:) = eigmax(:);
%     eigmins(i_Lelm,:) = eigmin(:);
%     eigmaxes_speckle(i_Lelm,:) = eigmax_speckle(:);
%     eigmins_speckle(i_Lelm,:) = eigmin_speckle(:);
%     eigmaxes_3pkt(i_Lelm,:) = eigmax_3pkt(:);
%     eigmins_3pkt(i_Lelm,:) = eigmin_3pkt(:);
    %%

    focus_ind_nonfocus = find(abs(depth_axis-20*1e-3)<1e-4);
    focus_ind_nonfocus = focus_ind_nonfocus(round(length(focus_ind_nonfocus)/2));
    
    d_orig_nonfocus = focus_ind_nonfocus:N:N*E;
    a_orig_nonfocus = N*az_inds:1:N*(az_inds+1);
    
    
    
    square_mid_ind_nonfocus = intersect(a_orig_nonfocus, d_orig_nonfocus);
    
    allinds_nonfocus = [];
    for i =1:25
        allinds_nonfocus = [allinds_nonfocus (square_mid_ind_nonfocus -(i-1) - length(depth_axis)*15):length(depth_axis):(square_mid_ind_nonfocus - (i-1) + length(depth_axis)*15)];
        allinds_nonfocus = [allinds_nonfocus (square_mid_ind_nonfocus +(i-1) - length(depth_axis)*15):length(depth_axis):(square_mid_ind_nonfocus + (i-1) + length(depth_axis)*15)];
    end
    allinds_nonfocus = sort(allinds_nonfocus);
    
    tmpstruct = mv_getCapon_speckle.StructTmp(:);
    
    nBins = 100;
    tmpstruct_ev_nonfocus = tmpstruct(allinds_nonfocus);%.R_eigvals
    
    
    hist_edges_nonfocus = linspace(0,2e-3, nBins+1);
    
    [tothist_nonfocus, hist_edges_nonfocus] = histcounts(tmpstruct_ev_nonfocus(1).R_eigvals_DL, nBins, 'BinEdges', hist_edges_nonfocus);
    for i = 2:length(tmpstruct_ev_nonfocus)
    %     if i ~= i_save
            tothist_nonfocus = tothist_nonfocus + histcounts(tmpstruct_ev_nonfocus(i).R_eigvals_DL, nBins, 'BinEdges', hist_edges_nonfocus);
    %     end
    end
        
    
    figure(1234+i_Lelm)
    subplot(211)
    h1 = histogram('BinCounts', tothist_nonfocus, 'BinEdges', hist_edges_nonfocus);
    xlabel('lambda')
    ylabel('Ant. lambda')
    title('Outside of focus')
    
    
    % %%
    
    focus_ind_focus = find(abs(depth_axis-40*1e-3)<1e-4);
    focus_ind_focus = focus_ind_focus(round(length(focus_ind_focus)/2));
    
    d_orig_focus = focus_ind_focus:N:N*E;
    a_orig_focus = N*az_inds:1:N*(az_inds+1);
    
    square_mid_ind_focus = intersect(a_orig_focus, d_orig_focus);
    
    allinds_focus = [];
    for i =1:25
        allinds_focus = [allinds_focus (square_mid_ind_focus -(i-1) - length(depth_axis)*15):length(depth_axis):(square_mid_ind_focus - (i-1) + length(depth_axis)*15)];
        allinds_focus = [allinds_focus (square_mid_ind_focus +(i-1) - length(depth_axis)*15):length(depth_axis):(square_mid_ind_focus + (i-1) + length(depth_axis)*15)];
    end
    allinds_focus = sort(allinds_focus);
    tmpstruct_ev_focus = tmpstruct(allinds_focus);%.R_eigvals
    
    hist_edges_focus = linspace(0,2e-4, nBins+1);
    
    [tothist_focus, hist_edges_focus] = histcounts(tmpstruct_ev_focus(1).R_eigvals_DL, nBins, 'BinEdges', hist_edges_focus);
    
    for i = 2:length(tmpstruct_ev_focus)
    %     if i ~= i_save
            tothist_focus = tothist_focus + histcounts(tmpstruct_ev_focus(i).R_eigvals_DL, nBins, 'BinEdges', hist_edges_focus);
    %     end
    end
        
    
    subplot(212)
    h2 = histogram('BinCounts', tothist_focus, 'BinEdges', hist_edges_focus);
    xlabel('lambda')
    ylabel('Ant. lambda')
    title('At focus')

    histcounts_tot_foc(i_Lelm, :) = tothist_focus;
    histcounts_tot_nonfoc(i_Lelm, :) = tothist_nonfocus;
    %%
    
end
%%
% DAStmp = b_data_DAS.data(:);
% tmp1 = depth_inds:length(depth_axis):length(DAStmp);
% tmp2 = (depth_inds-1):length(depth_axis):length(DAStmp);
% tmp3 = (depth_inds-2):length(depth_axis):length(DAStmp);
% tmp4 = (depth_inds+1):length(depth_axis):length(DAStmp);
% tmp5 = (depth_inds+2):length(depth_axis):length(DAStmp);
% for ind = 1:length(tmp1)
%     tmpresind = sort([tmp1(ind) tmp2(ind) tmp3(ind) tmp4(ind) tmp5(ind)]);
%     TotRes(1, ind) = max(abs(DAStmp(tmpresind)));
%     TotRes_3pkt(1, ind) = max(abs(b_data_DAS_3pkt.data(tmpresind)));
%     TotRes_speckle(1, ind) = max(abs(b_data_DAS_speckle.data(tmpresind)));
% end   
% LelmRes(1,:) = TotRes(1,:);
% 
% maxDas          = max(abs(DAStmp(:)));
% maxDas_3pkt     = max(abs(b_data_DAS_3pkt.data(:)));
% maxDas_speckle  = max(abs(b_data_DAS_speckle.data(:)));
% 
% %%
% % if ~exist(fullfile(basepath, 'Data', date)); mkdir(fullfile(basepath, 'Data', date)); end
% % save(fullfile(basepath, 'Data', date,'Lelm_analysis.mat'), 'TotRes', 'TotRes_Pow', 'LelmRes', 'LelmRes_Pow', 'maxDas', 'Lelms', 'CN_mid', 'CN_max', 'CN_left', 'CN_right', 'CN_mid_above', 'CN_mid_below', 'CN_rand_above', 'CN_rand_below', 'CN_Speckle_results')
% %% Result plots
% b_data_mv_getCapon.plot(figure(100), 'Capon on full scene')
% b_data_mv_getCapon_3pkt.plot(figure(101), 'Capon on only 3 points')
b_data_mv_getCapon_Speckle.plot(figure(102), 'Capon on speckle')


%% Resolution plots
% figure(110)
% subplot(311)
% plot(x_axis_dmm, db(abs(TotRes./maxDas)), 'LineWidth', 1.25);
% xlabel('x [mm]')
% ylabel('Amplitude [dB]')
% title('Full scene')
% legend('DAS', 'Capon')
% xlim([-15 15])
% ylim([-60 5])
% yticks(-66:12:6)
% grid on
% 
% subplot(312)
% plot(x_axis_dmm, db(abs(TotRes_3pkt./maxDas_3pkt)), 'LineWidth', 1.25);
% xlabel('x [mm]')
% ylabel('Amplitude [dB]')
% title('Scatterer scene')
% legend('DAS', 'Capon')
% xlim([-15 15])
% ylim([-60 5])
% yticks(-66:12:6)
% grid on
% 
% subplot(313)
% plot(x_axis_dmm, db(abs(TotRes_speckle./maxDas_speckle)), 'LineWidth', 1.25);
% xlabel('x [mm]')
% ylabel('Amplitude [dB]')
% title('Speckle scene')
% legend('DAS', 'Capon')
% xlim([-15 15])
% % ylim([-60 5])
% yticks(-60:6:6)
% grid on



%%



%% 
% idxes = B(setdiff(1:end,allinds_focus))
copy = uff.beamformed_data(b_data_mv_getCapon_Speckle);
copy.data(setdiff(1:end,[allinds_focus allinds_nonfocus])) = 0.01;
copy.plot(figure(1000))



%%
figure(123)
subplot(211)
plot(linspace(hist_edges_focus(1), hist_edges_focus(end),nBins), 1./histcounts_tot_foc)
title('Total histcounts inside focus')
subplot(212)
plot(linspace(hist_edges_nonfocus(1), hist_edges_nonfocus(end),nBins), 1./histcounts_tot_nonfoc)
title('Total histcounts outside focus')



