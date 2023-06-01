%% eigval_scripts.m
% Script to generate points at a distance d from each other.
clear;
close all;

pos = 20;%[20, 30, 40];
dist = [0, 2, 5, 10, 20];


for p_i = 1:length(pos)
    for d_i = 1:length(dist)
        d = dist(d_i);
        p = pos(p_i);
        
        if d == 0
            phantom_positions = [0, 0, p*1e-3];
        else
            phantom_positions = zeros(2, 3);
            
            phantom_positions(1,:) = [0.5*d*1e-3, 0, p*1e-3];
            phantom_positions(2,:) = [-0.5*d*1e-3, 0, p*1e-3];

        end
        
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
        regCoeff = 1/100;
        
        if ~exist('Lelm', 'var')
            Lelm = channel_data.probe.N*L_frac;
        end
        
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
        scan=uff.linear_scan('azimuth_axis',azimuth_axis,'depth_axis',depth_axis);
        
        
        script_mid_DAS_MLA
        
        
        %% ustb getcapon postprocess
        Lelm_set=0;
        calc_param = 0;
        script_post_getCapon
        b_data_mv_getCapon = mv_getCapon.go();
        
        %% Plot eigenvalues
        N = mv_getCapon.scan.N_depth_axis;
        E = mv_getCapon.scan.N_azimuth_axis;
        az = 0; % Rad
        az_inds = find(abs(azimuth_axis-az)<0.01);
        az = 0; % Degrees
        az_inds = find(abs(azimuth_axis-deg2rad(az))<0.01);
        az_inds = az_inds(round(length(az_inds)/2));

        
        fighandle1 = figure;
%         for i = 1:length(depth_axis)
%             for j = az_inds
%                 Rtmp = mv_getCapon.StructTmp(i,j).R;
%                 eigvals = eig(Rtmp);
%                 plot(depth_axis(i)*ones(length(eigvals),1)*1e3, eigvals, 'o')
%                 hold on
%             end
%         end
%         hold off
%         ylabel("Eigenvalues")
%         xlabel("Depth [mm]")
%         title("Eigenvalues along " + num2str(az) + " degrees")
%         
        max_length = max(max(cellfun('length',{mv_getCapon.StructTmp.R_eigvals})));
        Reig_axial = zeros(N, max_length);
        
        for i = 1:length(depth_axis)
            for j = az_inds
                len = cellfun('length',{mv_getCapon.StructTmp(i,j).R_eigvals});
                Reig_axial(i, 1:len) = sort(mv_getCapon.StructTmp(i,j).R_eigvals(:), 'descend');
                Reig_axial(i,(len+1):end) = -0.1;
            end
        end
        
        
%         figure(100)
        subplot(121)
        imagesc(Reig_axial.');
        ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
        xlabel("Depth [pixels]")
        colbar = colorbar();
        colbar.Title.String = "Eigenvalue";
        subplot(122)
        imagesc(db(Reig_axial).')
        colbar = colorbar();
        colbar.Title.String = "Eigenvalue [dB]";
        clim([-60 0])
        ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
        xlabel("Depth [pixels]")
        
        sgtitle("Eigenvalues sorted along " + num2str(az) + " degrees, NaN showed as value -0.1 [-20dB]")

    
%         savefig(fighandle1, ...
%             "eigvals_"+num2str(d)+"dist_"+num2str(p)+"depth", ...
%             path = "Eigenvalues\First test", overwrite = 1)
        %% Plot condition number
%         N = mv_getCapon.scan.N_depth_axis;
%         E = mv_getCapon.scan.N_azimuth_axis;
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
                [CN(i,j), CN_DL(i,j), Tse(i,j), Reigval(i,j), Reigval_DL(i,j), R_eigval] = parameter_calculations(R_tmp, R_DL_tmp, w_tmp, M_new_tmp);
            end
        end
        b_data_CondNr                =   uff.beamformed_data(b_data_mv_getCapon);
        b_data_CondNr_DL             =   uff.beamformed_data(b_data_mv_getCapon);
        b_data_CondNr.data(:,:)      =   CN(:);
        b_data_CondNr_DL.data(:,:)   =   CN_DL(:);
        fighandle2 = figure;
        b_data_CondNr.plot(fighandle2, ['Condition Number, d = ', num2str(d), 'depth = ',num2str(p)], [], 'none')
        fighandle3 = figure;
        b_data_CondNr_DL.plot(fighandle3, ['Condition Number, d = ', num2str(d), 'depth = ',num2str(p)], [], 'none')
 
%         savefig(fighandle2, ...
%             "CN_"+num2str(d)+"dist_"+num2str(p)+"depth", ...
%             path = "Condition Number\First test", overwrite = 1)
%         savefig(fighandle3, ...
%             "CNDL_"+num2str(d)+"dist_"+num2str(p)+"depth", ...
%             path = "Condition Number\First test", overwrite = 1)
    end
end

