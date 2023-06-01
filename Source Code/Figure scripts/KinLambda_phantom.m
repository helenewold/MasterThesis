%% Script explanation
% % Condition number output fra datasettet FieldII_noSpeckle4scats_5degs.uff
%
% Oppdatert: 21.11.22
close all
clear
%% Parameter creation
path = 'Datasets';
file_name = 'P4_FI_121444_45mm_focus';

filename = append(fullfile(path, file_name),'.uff');

save_fig = 0;

receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;
MLA = 1;
MLA_overlap = 1;

K = 0:5:20;
L_frac = 1/3;
regCoeff = 1/100;


%% Read data
script_ReadData

%% Midprocess - Delay and sum
script_mid_DAS_MLA

%% (optional) reshape vekk til 1 frame
data_cube = USTB_reshape(scan_MLA, b_data_MLA_delayed, channel_data, n_frame = 1);


for k = 1:length(K)
    K_in_lambda = K(k);
    %% getCapon postprocess structure
    Lelm = channel_data.probe.N*L_frac;
    
    
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
    fig1 = figure(10+k);
    b_data_mv_getCapon.plot(fig1, ['getCapon, Ksamps = ', num2str(mv_getCapon.K_samples)]);
    
    %% Condition number plotting
    b_data_CondNr                   =   uff.beamformed_data(b_data_mv_getCapon);
    b_data_CondNr_after             =   uff.beamformed_data(b_data_mv_getCapon);
    
    b_data_CondNr.data(:,:)         =   mv_getCapon.CN;
    b_data_CondNr_after.data(:,:)   =   mv_getCapon.CN_after;
    
    %%
    fig2 = figure(20+k);
    fig2.WindowState = 'maximized';
    sgtitle(['Condition numbers, Ksamps = ', num2str(mv_getCapon.K_samples)])
    
    subplot(121)
    plot_CondNr(b_data_CondNr, azimuth_axis, depth_axis, title = "Before diagonal loading")
    
    
    subplot(122)
    plot_CondNr(b_data_CondNr_after, azimuth_axis, depth_axis, title = "After diagonal loading")
%%
    fig3 = figure(30+k);
    fig3.WindowState = 'maximized';
    sgtitle(['Condition numbers, log10, Ksamps = ', num2str(mv_getCapon.K_samples)])
    
    subplot(121)
    plot_CondNr(b_data_CondNr, azimuth_axis, depth_axis, title = "Before diagonal loading", dB = 2)
    
    
    subplot(122)
    plot_CondNr(b_data_CondNr_after, azimuth_axis, depth_axis, title = "After diagonal loading", dB = 2)
%%
    if exist('save_fig','var') && save_fig
        savefig(fig1, "getCapon_" + num2str(mv_getCapon.K_samples) + 'Ksamps_'     + file_name,     path = "getCapon_Outputs")

        savefig(fig2, "CN_" + num2str(mv_getCapon.K_samples) + 'Ksamps_'     + file_name,     path = "Condition_Number")
        savefig(fig3, "CN_log10_" + num2str(mv_getCapon.K_samples) + 'Ksamps_'     + file_name,     path = "Condition_Number")
    end

end

