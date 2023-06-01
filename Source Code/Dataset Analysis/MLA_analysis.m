%% MLA_analysis.m
% Script to generate points at a distance d from each other.
clear;
close all;
%%
if ispc
    addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\downloads\Field_II_ver_3_30_windows.tar'));
    addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\Dokumenter\USTB'));
    addpath(genpath('\\hume.uio.no\student-u55\helenewo\MasterThesis'));
    path_ = ['..\Figures\Server\',date];
    path_fig = ['\FiguresFigFormat\Server\',date];
    basepath = '\\hume.uio.no\student-u55\helenewo\MasterThesis\MasterThesis';
elseif isunix
    addpath(genpath('/hom/dsb/field'));
    addpath(genpath('/uio/hume/student-u55/helenewo/pc/Dokumenter/USTB'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/datasets'));
    path_ = ['./Figures/Server/',date];
    path_fig = ['./FiguresFigFormat/Server/',date];
    basepath = '/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis';
end

%%

p = 30;%[20, 30, 40];
d = 2;%[0, 2, 5, 10, 20];
degs = 360 * (d / (2*pi*p));

MLAs = [1 3 5];
TotRes = zeros(length(MLAs)*2, 128*max(MLAs));
az_lens = zeros(length(MLAs),1);
x_ax_store = zeros(length(MLAs), 128*max(MLAs));
%%
for MLAind = 1:length(MLAs)
    MLA = MLAs(MLAind);
    phantom_positions = zeros(2, 3);
    
    phantom_positions(1,:) = [ p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
    phantom_positions(2,:) = [-p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
    h.save = 0;
    h.amp = "ones";
    
    channel_data = fieldii_generate_dataset("tmp", "tmp", phantom_positions,h);
    %% Beamforming section - Kun getCapon beamforming
    % Parameter creation
    receive_window = uff.window.boxcar;
    f_number = 1.7;
    pw_margin = 5e-3;
    spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
    transmit_window = uff.window.scanline;
%     MLA = 1;
    MLA_overlap = 1;
    
    K_in_lambda = 2;
    % Lelm = channel_data.probe.N/3;
    L_frac = 1/3;
    regCoeff = 1/100;
    
    Lelm = channel_data.probe.N*L_frac;

    %% Read data seksjon
    % Setter altså variabler slik som i read data scriptet, men uten å lese fra
    % fil.
            
    channel_data.N_frames = 1;
            
    % Definerer scan
    azimuth_axis=zeros(channel_data.N_waves,1);
    depth_axis = linspace(1e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';
    
    for n=1:channel_data.N_waves
        azimuth_axis(n)=channel_data.sequence(n).source.azimuth;
    end
    script_mid_DAS_MLA
    script_mid_DAS_bothDims
    b_data_DAS = mid_DAS.go();

    x_axis_dmm = p*azimuth_MLA;

    %% ustb getcapon postprocess
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
    
    depth_ind = find(abs(depth_axis-p*1e-3)<1e-4);
    depth_ind = depth_ind(round(length(depth_ind)/2));


%             axises = zeros(N*E,1);
%             axises(depth_inds:N:end) = 1;           % Depth
%             axises(N*az_inds:1:N*(az_inds+1)) = 1;  % Azimuth
%             inds = find(axises == 1);
%             b_data_mv_getCapon.data(inds) = 10;
% 
% 
%     fig_getCapon = figure();
%     b_data_mv_getCapon.plot(fig_getCapon, ['getCapon output, MLA = ', num2str(MLA)])
%     xlim([-10 10])
%     ylim([25 35])
%     save_fig(fig_getCapon, ...
%         ['result_',num2str(MLA),'MLA'], ...
%         path = fullfile(path_, 'MLA'))

    tmp1 = depth_ind:length(depth_axis):length(b_data_mv_getCapon.data);
    tmp2 = (depth_ind-1):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp3 = (depth_ind-2):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp4 = (depth_ind+1):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp5 = (depth_ind+2):length(depth_axis):length(b_data_mv_getCapon.data);
    
    
    for ind = 1:length(tmp1)
        tmpresind = sort([tmp1(ind) tmp2(ind) tmp3(ind) tmp4(ind) tmp5(ind)]);
%         tmpRES(ind) = max(abs(b_data_mv_getCapon.data(tmpresind)));

        TotRes(MLAind, ind) = max(abs(b_data_mv_getCapon.data(tmpresind)));
        TotRes(MLAind+3, ind) = max(abs(b_data_DAS.data(tmpresind)));
    end
    tmp = b_data_mv_getCapon;
    tmp.data(tmp1) = 1;
    tmp.plot(figure)
    figure()
    plot(db(abs(TotRes(MLAind,:))))
    hold on
    plot(db(abs(TotRes(MLAind+3,:))))
    hold off
    az_lens(MLAind) = length(x_axis_dmm);
    x_ax_store(MLAind,1:length(x_axis_dmm)) = x_axis_dmm;
end
%% fig_MLA = figure(100);
% drawnow
% plot(x_axis_dmm, db(abs(tmpRES)), 'DisplayName', ['MLA = ', num2str(MLA)]);
% drawnow
% plot(x_axis_dmm, -6*ones(length(x_axis_dmm),1), '-r', 'DisplayName','-6dB')
% hold off
% title("Resolution through 30mm")
% xlabel("Width [mm]")
% ylabel("Amplitude [dB]")
% ylim([-60 5])
% xlim([x_axis_dmm(1) x_axis_dmm(end)])
% legend('Location', 'southeast')

% save_fig(fig_resolution, ...
%     "resolution_DL", ...
%     path = fullfile(date, 'Resolution'))
maxDAS = max(max(abs(TotRes(4:6, :))));
fig_MLA = figure(100);
fig_MLA.Position = [200 100 1000 400];
for MLAi = 1:length(MLAs)
    drawnow
    plot(x_ax_store(MLAi, 1:az_lens(MLAi)), db(abs(TotRes(MLAi, 1:az_lens(MLAi)))./maxDAS), 'LineWidth',1.25,'DisplayName', ['MV, MLA = ', num2str(MLAs(MLAi))]);
    hold on
    plot(x_ax_store(MLAi, 1:az_lens(MLAi)), db(abs(TotRes(MLAi+3, 1:az_lens(MLAi)))./maxDAS), '--','LineWidth',1.25, 'DisplayName', ['DAS, MLA = ', num2str(MLAs(MLAi))]);
    hold on
end
hold off
legend('Location', 'northeast')
title("Lateral line through 30[mm]")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 3])
% xlim([-8 8])
yticks(-60:6:6)
grid on
save_fig(fig_MLA, ...
    "resolution_MLA", ...
    path = fullfile(path_,'MLA'), fig_path=fullfile(path_, 'MLA'), base_path = "")