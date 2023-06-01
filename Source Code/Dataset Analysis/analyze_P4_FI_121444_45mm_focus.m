close all
clear

%% Parameter creation
path = 'Datasets';
file_name = 'P4_FI_121444_45mm_focus';

filename = append(fullfile(path, file_name),'.uff');

save_fig = 0;
ow = 1;

plt = 0;


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

%% Midprocess - DAS both dimensions
% b_data_DAS
script_mid_DAS_bothDims


%% (optional) reshape vekk til 1 frame
data_cube = USTB_reshape(scan_MLA, b_data_MLA_delayed, channel_data, n_frame = 1);

Lelm = channel_data.probe.N*L_frac;

%% USTB Capon
% b_data_mv_MLA
script_post_capon_minvar

%% getCapon postprocess
script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();

%% Tse first run part
b_data_Tse             =   uff.beamformed_data(b_data_mv_getCapon);
b_data_Tse.data(:,:)   =   mv_getCapon.Tse;


%% Condition number of this
b_data_CondNr                   =   uff.beamformed_data(b_data_mv_getCapon);
b_data_CondNr_after             =   uff.beamformed_data(b_data_mv_getCapon);

b_data_CondNr.data(:,:)         =   mv_getCapon.CN;
b_data_CondNr_after.data(:,:)   =   mv_getCapon.CN_after;


%% Diff image
b_data_diff = uff.beamformed_data(b_data_mv_MLA);

% data_diff = abs(b_data_mv_MLA.data(:) - post_getCapon.data(:));
b_data_diff.data(:,1,1,1)  = abs(abs(b_data_mv_MLA.data./max(b_data_mv_MLA.data(:))) ...
                        -abs(b_data_mv_getCapon.data./max(b_data_mv_getCapon.data(:))));



%% Histogram matching of this
[img_out, b_data_out] = tools.histogram_match(b_data_DAS,b_data_mv_getCapon);
    


%% change L
Lfrac = [3/4, 2/3, 1/2, 1/3, 1/4];
Lfrac = 1/3;
K_in_lambda = 0.5;

for i = 1:length(Lfrac)
    L_frac = Lfrac(i);
    Lelm = channel_data.probe.N*L_frac;

    script_post_getCapon
    b_data_Lfractest = mv_getCapon.go();
    
    disp(['L = ', num2str(L_frac), ', K = ', num2str(mv_getCapon.K_samples)])
    
    if plt
        fig_getCapon_Lfrac = figure(100+i);
        b_data_Lfractest.plot(fig_getCapon_Lfrac, ['getCapon, L = M*', num2str(L_frac)]);
    end

    b_data_CondNr_Lfrac                   =   uff.beamformed_data(b_data_mv_getCapon);
    b_data_CondNr_after_Lfrac             =   uff.beamformed_data(b_data_mv_getCapon);
    
    b_data_CondNr_Lfrac.data(:,:)         =   mv_getCapon.CN;
    b_data_CondNr_after_Lfrac.data(:,:)   =   mv_getCapon.CN_after;
    
    
%     fig_CondNr_Lfrac = figure(150);
%     fig_CondNr_Lfrac.WindowState = 'maximized';
%     sgtitle(['Condition numbers, L = M*' , num2str(L_frac)])
%     
%     subplot(121)
%     plot_CondNr(b_data_CondNr_Lfrac, azimuth_axis, depth_axis, title = "Before diagonal loading")
%     subplot(122)
%     plot_CondNr(b_data_CondNr_after_Lfrac, azimuth_axis, depth_axis, title = "After diagonal loading")
    
    fig_CondNrLog10_Lfrac = figure(200);
%     fig_CondNrLog10_Lfrac.WindowState = 'maximized';
%     sgtitle(['Condition numbers, log10, L = M*' , num2str(L_frac)])
%     subplot(121)
    plot_CondNr(b_data_CondNr_Lfrac, azimuth_axis, depth_axis, title = ['Condition numbers, log10, L = M*' , num2str(L_frac), "Before"], dB = 2)
%     subplot(122)
%     plot_CondNr(b_data_CondNr_after_Lfrac, azimuth_axis, depth_axis, title = "After diagonal loading", dB = 2)


%             if exist('save_fig','var') && save_fig
%                 savefig(fig_getCapon_Lfrac, "getCapon_" + num2str(L_frac) + 'M_' + file_name,     path = fullfile("fig_" + file_name,"getCapon_Outputs"), overwrite = ow)
%                 savefig(fig_CondNr_Lfrac, "CN_" + num2str(L_frac) + 'M_'         + file_name,     path = fullfile("fig_" + file_name,"Condition_Number"), overwrite = ow)
%                 savefig(fig_CondNrLog10_Lfrac, "CN_log10_" + num2str(L_frac) + 'M_'   + file_name,     path = fullfile("fig_" + file_name,"Condition_Number"), overwrite = ow)
%             end
%     end

end

%% change K
% if section_Ksamps
L_frac = 1/3;
Lelm = channel_data.probe.N*L_frac;

K = 0.5;        %0:0.5:5;

for k = 1:length(K)
    K_in_lambda = K(k);


    % Postprocess - getCapon
    script_post_getCapon
    b_data_Ksamps = mv_getCapon.go();

%     disp(['L = ', num2str(L_frac), ', K = ', num2str(mv_getCapon.K_samples)])

    % Plot result
    if plt 
        fig_getCapon_Ksamps = figure(500+k);
        b_data_Ksamps.plot(fig_getCapon_Ksamps, ['getCapon, Ksamps = ', num2str(mv_getCapon.K_samples)]);
    end
    % Condition number plotting
    b_data_CondNr_Ksamps                   =   uff.beamformed_data(b_data_mv_getCapon);
    b_data_CondNr_after_Ksamps             =   uff.beamformed_data(b_data_mv_getCapon);
    
    b_data_CondNr_Lfrac.data(:,:)         =   mv_getCapon.CN;
    b_data_CondNr_after_Lfrac.data(:,:)   =   mv_getCapon.CN_after;
    
    %
%     if plt
%     fig_CondNr_Ksamps.("K" + mv_getCapon.K_samples) = figure(550+k);
%     fig_CondNr_Ksamps.("K" + mv_getCapon.K_samples).WindowState = 'maximized';
%     sgtitle(['Condition numbers, Ksamps = ', num2str(mv_getCapon.K_samples)])
%     
%     subplot(121)
%     plot_CondNr(b_data_CondNr_Ksamps, azimuth_axis, depth_axis, title = "Before diagonal loading")
%     
%     
%     subplot(122)
%     plot_CondNr(b_data_CondNr_after_Ksamps, azimuth_axis, depth_axis, title = "After diagonal loading")

    fig_CondNrLog10_Ksamps.("K" + mv_getCapon.K_samples) = figure(600+k);
%     fig_CondNrLog10_Ksamps.("K" + mv_getCapon.K_samples).WindowState = 'maximized';
%     sgtitle(['Condition numbers, log10, Ksamps = ', num2str(mv_getCapon.K_samples)])
    
%     subplot(121)
    plot_CondNr(b_data_CondNr_Ksamps, azimuth_axis, depth_axis, title = ['Condition numbers, log10, Ksamps = ', num2str(mv_getCapon.K_samples), "Before"], dB = 2)
    
    
%     subplot(122)
%     plot_CondNr(b_data_CondNr_after_Ksamps, azimuth_axis, depth_axis, title = "After diagonal loading", dB = 2)
%     
%         if exist('save_fig','var') && save_fig
%             savefig(fig_getCapon_Ksamps,    "getCapon_Ksamps_" + num2str(mv_getCapon.K_samples) + file_name, path = fullfile("fig_" + file_name, "getCapon_Outputs"), overwrite = ow)
%             savefig(fig_CondNr_Ksamps,      "CN_Ksamps_"       + num2str(mv_getCapon.K_samples) + file_name, path = fullfile("fig_" + file_name,"Condition_Number"), overwrite = ow)
%             savefig(fig_CondNrLog10_Ksamps, "CN_log10_Ksamps_" + num2str(mv_getCapon.K_samples) + file_name, path = fullfile("fig_" + file_name,"Condition_Number"), overwrite = ow)
%         end
end
% end

%% Plot and save
if plt
    fig1 = figure(1);
    b_data_DAS.plot(fig1, 'DAS');
    fig2 = figure(2);
    b_data_mv_MLA.plot(fig2, 'MV Capon');
    fig3 = figure(3);
    b_data_mv_getCapon.plot(fig3, 'getCapon');

    figTse = figure(12345);
    b_data_Tse.plot(figTse, "Tse", [], 'none');

    fig4 = figure(4);
    fig4.WindowState = 'maximized';
    sgtitle('Condition numbers')
    subplot(121)
    plot_CondNr(b_data_CondNr, azimuth_axis, depth_axis, title = "Before diagonal loading", dB = 1)
    subplot(122)
    plot_CondNr(b_data_CondNr_after, azimuth_axis, depth_axis, title = "After diagonal loading", dB = 1)

    fig5 = figure(5);
    fig5.WindowState = 'maximized';
    sgtitle('Condition numbers - log10(data)')
    subplot(121)
    plot_CondNr(b_data_CondNr, azimuth_axis, depth_axis, title = "Before diagonal loading", dB=2)
    subplot(122)
    plot_CondNr(b_data_CondNr_after, azimuth_axis, depth_axis, title = "After diagonal loading", dB=2)

    fig6 = figure(6);
    b_data_diff.plot(fig6,'Difference between getCapon and capon\_minvar');

    fig7 = figure(7);
    b_data_out.plot(fig7, "Histogram matching", 60, 'none');
    caxis([-60 0])
end


if exist('plt','var') && plt && exist('save_fig','var') && save_fig
    savefig(fig1, "DAS_"        + file_name,     path = "fig_" + file_name, overwrite = ow)
    savefig(fig2, "MVCapon_"    + file_name,     path = "fig_" + file_name, overwrite = ow)
    savefig(fig3, "getCapon_"   + file_name,     path = "fig_" + file_name, overwrite = ow)
    savefig(fig4, "CN_"         + file_name,     path = "fig_" + file_name, overwrite = ow)
    savefig(fig5, "CNlog10_"    + file_name,     path = "fig_" + file_name, overwrite = ow)
    savefig(fig6, "diff_"       + file_name,     path = "fig_" + file_name, overwrite = ow)
    savefig(fig7, "HM_"         + file_name,     path = "fig_" + file_name, overwrite = ow)
    savefig(figTse, "Tse_"    + file_name,     path = "fig_" + file_name, overwrite = ow)
end
