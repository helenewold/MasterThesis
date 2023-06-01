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
else
    basepath = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode';
end
cd(fullfile(basepath, 'Source Code'))
% if ~exist('speckle_value')
speckle_value=0.5;
% end
data_path = fullfile(basepath, 'Data', 'Fullrun', ['SpeckleLevel', num2str(speckle_value)], 'ddist');
if ~exist(data_path); mkdir(data_path); end
evpath = fullfile(basepath, 'Figures', 'Fullrun', ['SpeckleLevel', num2str(speckle_value)], 'ddistEV');
if ~exist(evpath); mkdir(evpath); end
path_ = fullfile(basepath, 'Figures', 'Fullrun',['SpeckleLevel', num2str(speckle_value)], 'Ddist');
if ~exist(path_); mkdir(path_); end

%% Load Data seksjon
% last inn speckle datasett
filename = fullfile('datasets', 'speckle_simulation_v1.uff');

receive_window = uff.window.boxcar;
f_number = 1.7;
pw_margin = 5e-3;
spherical_transmit_delay_model_ = spherical_transmit_delay_model.hybrid;
transmit_window = uff.window.scanline;



MLA = 3;
MLA_overlap = 1;

% Capon values (Kan endres underveis i scriptet)
K_in_lambda = 2;
% L_frac = 1/3;
regCoeff = 1/100;
%% Read data
script_ReadData
channel_data_speckle = channel_data;

%% DAS transmit dimension
% b_data_MLA_delayed
script_mid_DAS_MLA
speckle_b_data = b_data_MLA_delayed;
%% Kombiner datasett


p = 30;%[20, 30, 40];
dist = [0, 0.5, 1, 1.5, 2, 5];

degs_separated = 360 * (dist / (2*pi*p));

L_frac = 1/3;

DistResMVAmp = zeros(length(dist), 128*MLA);
DistResMVPow = zeros(length(dist), 128*MLA);
DistResDAS = zeros(length(dist), 128*MLA);
maxDasVals = zeros(length(dist), 1);
ii=1;

for degs_i = 1:length(degs_separated)
    disp(['dist no ', num2str(degs_i), 'of ', num2str(length(degs_separated))])
    degs = degs_separated(degs_i);
    d = dist(degs_i);

    phantom_positions(1,:) = [ p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
    phantom_positions(2,:) = [-p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
        
    h.save = 0;
    h.amp = "ones";
    
    channel_data = fieldii_generate_dataset("tmp", "tmp", phantom_positions,h);

     %% Beamforming section - Kun getCapon beamforming
    % Parameter creation

    Lelm = channel_data.probe.N*L_frac;
    %% Read data seksjon
    % Setter altså variabler slik som i read data scriptet, men uten å lese fra
    % fil.
    
    channel_data.N_frames = 1;
    

    script_mid_DAS_MLA
    pkt_b_data = b_data_MLA_delayed;

    b_data_MLA_delayed.data(:,:) = pkt_b_data.data(:,:) + speckle_value*speckle_b_data.data(:,:);

    cc = postprocess.coherent_compounding();
    cc.input = b_data_MLA_delayed;
    b_data_DAS = cc.go();
    
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
    
    depth_inds = find(abs(depth_axis-p*1e-3)<1e-4);
    depth_inds = depth_inds(round(length(depth_inds)/2));

    fighandle4 = figure();
    b_data_mv_getCapon.plot(fighandle4, ['getCapon Output, d = ', num2str(d)])


    b_data_mv_powCapon = uff.beamformed_data(b_data_MLA_delayed);
    b_data_mv_powCapon.data = mv_getCapon.imPower;
    fig_powCap = figure();
    b_data_mv_powCapon.plot(fig_powCap, ['Power Capon output, d = ', num2str(d)])

    save_fig(fighandle4, ...
        "result_"+num2str(d)+"dist_"+num2str(p)+"depth", ...
        path = path_, fig_path=path_, base_path='')
    save_fig(fig_powCap, ...
        "result_powCap_"+num2str(d)+"dist_"+num2str(p)+"depth", ...
        path = path_, fig_path=path_, base_path='')

    %%
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
%%
    fighandle1 = figure; 

    drawnow
    imagesc(depth_axis*1e3, 1:max_length_az, db(Reig_axial).');
    ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
    xlabel("Depth [mm]")
    title("Eigenvalues of R through " + num2str(az) + " degrees, d="+num2str(d))
    colbar = colorbar();
    colbar.Title.String = "Eigenvalue [dB]";
    xlim([20 40])
    fighandle2 = figure; 
    drawnow
    imagesc(x_axis_dmm, 1:max_length_depth, db(Reig_depth).')
    colbar = colorbar();
    colbar.Title.String = "Eigenvalue [dB]";
    title("Eigenvalues of R through 30 mm, d="+num2str(d))
    ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
    xlabel("x [mm]")

    savefig(fighandle1, fullfile(evpath, ['EV_',num2str(d),'dist_azim_speckle.fig']))
    savefig(fighandle2, fullfile(evpath, ['EV_',num2str(d),'dist_depth_speckle.fig']))

    %% Gather resolution
    tmp1 = depth_inds:length(depth_axis):length(b_data_mv_getCapon.data);
    tmp2 = (depth_inds-1):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp3 = (depth_inds-2):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp4 = (depth_inds+1):length(depth_axis):length(b_data_mv_getCapon.data);
    tmp5 = (depth_inds+2):length(depth_axis):length(b_data_mv_getCapon.data);
    for ind = 1:length(tmp1)
        tmpresind = sort([tmp1(ind) tmp2(ind) tmp3(ind) tmp4(ind) tmp5(ind)]);
        DistResMVAmp(ii, ind) = max(abs(b_data_mv_getCapon.data(tmpresind)));
        DistResMVPow(ii, ind) = max(abs(mv_getCapon.imPower(tmpresind)));
        DistResDAS(ii, ind) = max(abs(b_data_DAS.data(tmpresind)));
    end

    maxDasVals(ii) = max(b_data_DAS.data(:,:));
    ii=ii+1;
end
 

save(fullfile(data_path, 'ddist_speckle.mat'), 'maxDasVals', 'DistResMVAmp', 'DistResMVPow','DistResDAS', "depth_axis", 'azimuth_MLA', 'x_axis_dmm')

%%
maxDas = maxDasVals(end);

DistRes(ii,:) = b_data_DAS.data(depth_inds:N:end);

L = {['d=',num2str(dist(1)), '[mm]'],...
    ['d=',num2str(dist(2)), '[mm]'],...
    ['d=',num2str(dist(3)), '[mm]'],...
    ['d=',num2str(dist(4)), '[mm]'],...
    ['d=',num2str(dist(5)), '[mm]'],...
    ['d=',num2str(dist(6)), '[mm]']
    };
lt = {'-', ':', ':', '--', '--','-'};

figure_res_MVAmp = figure(1234);
drawnow;
figure_res_MVAmp.Position = [200 100 1000 400];
for iii = 1:length(dist)
    plot(x_axis_dmm, db(abs(DistResMVAmp(iii,:))./maxDas), 'LineWidth',1.25, 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii)))
%     plot(x_axis_dmm, db(abs(DistResMVAmp(iii,:))), 'LineWidth',1.25, 'LineStyle', lt(iii),'DisplayName', string(L(iii)))
    hold on
end
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm] - Amplitude Capon")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 6])
xlim([-6 6])
yticks(-60:6:60)
grid on

%%
figure_res_MVPow = figure(1235);
drawnow;
figure_res_MVPow.Position = [200 100 1000 400];
for iii = 1:length(dist)
    plot(x_axis_dmm, db(abs(DistResMVPow(iii,:))), 'LineWidth',1.25, 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii)))
    hold on
end
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm] - Power Capon")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([max(max(db(abs(DistResMVPow))))-60 max(max(db(abs(DistResMVPow))))+6])
xlim([-6 6])
yticks(-84:6:60)
grid on

%%

figure_res_DAS = figure(1236);
drawnow;
figure_res_DAS.Position = [200 100 1000 400];
for iii = 1:length(dist)
    plot(x_axis_dmm, db(abs(DistResDAS(iii,:)./maxDas)), 'LineWidth',1.25, 'LineStyle', string(lt(iii)),'DisplayName', string(L(iii)))
    hold on
end
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm] - DAS")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 6])
xlim([-6 6])
yticks(-60:6:6)
grid on

%% Save
save_fig(figure_res_MVAmp, 'resolution_AmpCapon', path = path_, fig_path =path_, base_path='')
save_fig(figure_res_MVPow, 'resolution_PowCapon', path = path_, fig_path =path_, base_path='')
save_fig(figure_res_DAS,   'resolution_DAS', path = path_, fig_path =path_, base_path='')

