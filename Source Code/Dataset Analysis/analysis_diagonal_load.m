%% diagonal load.m
% Script to generate points at a distance d from each other.
clear;
close all;
%%
if isunix
    addpath(genpath('/hom/dsb/field'));
    addpath(genpath('/uio/hume/student-u55/helenewo/pc/Dokumenter/USTB'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/datasets'));
    path_ = ['./Figures/Server/',date];
    path_fig = ['./FiguresFigFormat/Server/',date];
    basepath = '/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis';
end
%%

pos = 30;%[20, 30, 40];
dist = 2;%[0, 2, 5, 10, 20];
degs_separated = 360 * (dist / (2*pi*pos));

% DL = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 1.5, 2, 2.5, 5:5:100]/100;
DL = 1/100;
CN_max = zeros(length(DL), 1);
%%
% DL = [0,  1/100, 1];
for p_i = 1:length(pos)
    for degs_i = 1:length(degs_separated)
        degs = 0.5*degs_separated(degs_i);
        d = dist(degs_i);
        p = pos(p_i);
%         if d == 0
%             phantom_positions = [0, 0, p*1e-3];
%         else
        phantom_positions = zeros(2, 3);
        
        phantom_positions(1,:) = [ p*1e-3*sind( degs ), 0, p*1e-3*cosd( degs )];
        phantom_positions(2,:) = [-p*1e-3*sind( degs ), 0, p*1e-3*cosd( degs )];

%         end
            
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
        MLA = 1;
        MLA_overlap = 1;
        
        K_in_lambda = 2;
        % Lelm = channel_data.probe.N/3;
        L_frac = 1/3;
%             regCoeff = 1/100;
        
        if ~exist('Lelm', 'var')
            Lelm = channel_data.probe.N*L_frac;
        end
        for DL_i = 1:length(DL)
            tic
            regCoeff = DL(DL_i);
            
            %% Read data seksjon
            % Setter altså variabler slik som i read data scriptet, men uten å lese fra
            % fil.
            
            channel_data.N_frames = 1;
            
            % Definerer scan
            azimuth_axis=zeros(channel_data.N_waves,1);
            depth_axis = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';
            
            for n=1:channel_data.N_waves
                azimuth_axis(n)=channel_data.sequence(n).source.azimuth;
            end
            script_mid_DAS_MLA
            
            x_axis_dmm = p*(azimuth_axis);

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


%             axises = zeros(N*E,1);
%             axises(depth_inds:N:end) = 1;           % Depth
%             axises(N*az_inds:1:N*(az_inds+1)) = 1;  % Azimuth
%             inds = find(axises == 1);
%             b_data_mv_getCapon.data(inds) = 10;
% 
% 
            fig_getCapon = figure();
            b_data_mv_getCapon.plot(fig_getCapon, ['getCapon Output, DL = ', num2str(regCoeff*100), '%'])
            xlim([-10 10])
            ylim([25 35])
%             save_fig(fig_getCapon, ...
%                 ['result_',num2str(regCoeff),'DL'], ...
%                 path = fullfile(date, 'Result'))

            max_length_az = max(cellfun('length',{mv_getCapon.StructTmp(:,az_inds).R_eigvals_DL}));
            max_length_depth = max(cellfun('length',{mv_getCapon.StructTmp(depth_inds,:).R_eigvals_DL}));
            Reig_axial = zeros(N, max_length_az);
            Reig_depth = zeros(E, max_length_depth);
            
            for i = 1:N
                len = cellfun('length',{mv_getCapon.StructTmp(i,az_inds).R_eigvals_DL});
                Reig_axial(i, 1:len) = sort(mv_getCapon.StructTmp(i,az_inds).R_eigvals_DL(:), 'descend');
                Reig_axial(i,(len+1):end) = nan;
            end
            for i = 1:E
                len = cellfun('length',{mv_getCapon.StructTmp(depth_inds,i).R_eigvals_DL});
                Reig_depth(i, 1:len) = sort(mv_getCapon.StructTmp(depth_inds,i).R_eigvals_DL(:), 'descend');
                Reig_depth(i,(len+1):end) = nan;
            end

            max_ax = max(max(Reig_axial));
            Reig_axial(isnan(Reig_axial)) = -0.1;

            max_dp = max(max(Reig_depth));
            Reig_depth(isnan(Reig_depth)) = -0.1;
            
%             fig_EV_azim = figure(); 
%             drawnow
%             imagesc(depth_axis*1e3, 1:max_length_az, db(Reig_axial).');
%             ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
%             xlabel("Depth [mm]")
%             title("Eigenvalues of R through " + num2str(az) + " degrees, DL="+num2str(regCoeff*100)+ '%')
%             colbar = colorbar();
%             colbar.Title.String = "Eigenvalue [dB]";
%             clim([db(max_ax)-60 db(max_ax)])
%             xlim([20 40])
% 
%             fig_EV_depth = figure; 
%             drawnow
%             imagesc(x_axis_dmm, 1:max_length_depth, db(Reig_depth).')
%             colbar = colorbar();
%             colbar.Title.String = "Eigenvalue [dB]";
%             clim([db(max_dp)-60 db(max_dp)])
%             title("Eigenvalues of R through 30 mm, DL="+num2str(regCoeff*100)+ '%')
%             ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
%             xlabel("Width [mm]")
% 
%             save_fig(fig_EV_azim, ...
%                 "EV_az_"+num2str(regCoeff)+"DL", ...
%                 path = fullfile(date, 'Eigenvalues'))
%             save_fig(fig_EV_depth, ...
%                 "EV_dp_"+num2str(regCoeff)+"DL", ...
%                 path = fullfile(date, 'Eigenvalues'))
            %% Plot condition number
            CN = zeros(N,E);
            CN_DL = zeros(N,E);
            Tse = zeros(N,E);
            Reigval = zeros(N,E);
            Reigval_DL = zeros(N,E);
            for i = 1:N
                for j = 1:E
                    R_tmp = mv_getCapon.StructTmp(i,j).R;
                    R_DL_tmp = mv_getCapon.StructTmp(i,j).R_DL;
                    w_tmp = mv_getCapon.StructTmp(i,j).w;
                    M_new_tmp = mv_getCapon.StructTmp(i,j).M;
                    [CN(i,j), CN_DL(i,j), Tse(i,j), Reigval(i,j), Reigval_DL(i,j)] = parameter_calculations(R_tmp, R_DL_tmp, w_tmp, M_new_tmp);
                end
            end
            b_data_CondNr_DL             =   uff.beamformed_data(b_data_mv_getCapon);
            b_data_CondNr_DL.data(:,:)   =   CN_DL(:);
%             fig_CN = figure;
%             b_data_CondNr_DL.plot(fig_CN, ['Condition Number, diagonal load ', num2str(regCoeff*100), '%'], [], 'none')
%      
%             save_fig(fig_CN, ...
%                 "CNDL_"+num2str(regCoeff)+"DL", ...
%                 path = fullfile(date, 'Condition Number'))
            CN_notnan = b_data_CondNr_DL.data(~isinf(b_data_CondNr_DL.data(~isnan(b_data_CondNr_DL.data))));
            CN_max(DL_i) = max(max(CN_notnan));
            disp(['DL no. ', num2str(DL_i), ' of ', num2str(length(DL)), ' done calculating'])
            CN_max(DL_i)
            fig_resolution = figure(101);
            fig_resolution.Position = [100 100 1000 400];
            drawnow
            plot(x_axis_dmm, db(abs(b_data_mv_getCapon.data(depth_inds:N:end)./max(b_data_mv_getCapon.data(depth_inds:N:end)))), 'DisplayName', ['DL=',num2str(regCoeff*100), '%'])
%             hold on
%             toc
        end
    end
end
% drawnow
% plot(x_axis_dmm, -6*ones(length(x_axis_dmm),1), '-r', 'DisplayName','-6dB')
% hold off
title("Resolution through 30mm")
xlabel("Width [mm]")
ylabel("Amplitude [dB]")
ylim([-10 5])
% xlim([x_axis_dmm(1) x_axis_dmm(end)])
% legend('Location', 'southeast')

% save_fig(fig_resolution, ...
%     "resolution_DL", ...
%     path = fullfile(date, 'Resolution'))
%% Plot max value of CN for all DLs
% figure_CNDL = figure(123);
% drawnow
% plot(DL, CN_max, '-o')
% xlabel("Diagonal load %")
% ylabel("max(Condition Number)")
% title("Maximum condition number for each diagonal load added to R")

% save_fig(figure_CNDL, "maxCNDL")

