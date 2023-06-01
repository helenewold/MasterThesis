%% eigval_scripts.m
% Script to generate points at a distance d from each other.
clear;
close all;

% kommentert ut fordi i tilfelle man ønsker å se på forskjellige
% plasseringer

pos = 30;%[20, 30, 40];
dist = [0, 2, 5, 10, 20];

Lfrac = 1/3;%[1/6, 1/4, 1/3, 1/2, 2/3, 3/4];
% fig_condnr_depth = figure(1000);
% fig_condnr_az = figure(1001);

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
            phantom_positions(i,:) = [pos(i)/1000*sind( degs_separated ), 0, pos(i)/1000*cosd( degs_separated )];
            phantom_positions(n_scats+i,:) = [-pos(i)/1000*sind( degs_separated ), 0, pos(i)/1000*cosd( degs_separated )];
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
        
        % Lelm = channel_data.probe.N/3;
        regCoeff = 1/100;
        K_in_lambda = 2;
        
        % Definerer scan
        azimuth_axis=zeros(channel_data.N_waves,1);
%         depth_axis = linspace(0e-3, 58e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';
        depth_axis = linspace(p*1e-3-10e-3, p*1e-3+10e-3, 512).'; %    z_axis=linspace(1e-3,55e-3,512).';

        for n=1:channel_data.N_waves
            azimuth_axis(n)=channel_data.sequence(n).source.azimuth;
        end  
            
        
        for Li = 1:length(Lfrac)
            L_frac = Lfrac(Li);
%             if ~exist('Lelm', 'var')
            Lelm = channel_data.probe.N*L_frac;
%             end
            %% Read data seksjon
            % Setter altså variabler slik som i read data scriptet, men uten å lese fra
            % fil.
            
            channel_data.N_frames = 1;
            

            script_mid_DAS_MLA
            
            %% ustb getcapon postprocess
            Lelm_set=0;
            calc_param = 0;
            script_post_getCapon
            b_data_mv_getCapon = mv_getCapon.go();
            %%
%             az = deg2rad(0); % deg
%             az_inds = find(abs(azimuth_axis-az)<0.01);
%             az_inds = az_inds(round(length(az_inds)/2));
%             depth = 30e-3;
%             depth_inds = find(abs(depth_axis-depth)<0.01);
%             depth_inds = depth_inds(round(length(depth_axis)/2));
 

            %% Plot eigenvalues
            K_samp = mv_getCapon.K_samples;
            N = mv_getCapon.scan.N_depth_axis;
            E = mv_getCapon.scan.N_azimuth_axis;
            az = deg2rad(0); % Rad
            az_inds = find(abs(azimuth_axis-az)<0.01);
            az_inds = az_inds(round(length(az_inds)/2));
            
            depth_inds = find(abs(depth_axis-30e-3)<1e-4);
            depth_inds = depth_inds(round(length(depth_inds)/2));
           
            axises = zeros(N*E,1);
            axises(depth_inds:N:end) = 1;           % Depth
            axises(N*az_inds:1:N*(az_inds+1)) = 1;  % Azimuth
            inds = find(axises == 1);
%             b_data_mv_getCapon.data(inds) = 10;
            fighandle4 = figure();
            b_data_mv_getCapon.plot(fighandle4, ['getCapon Output, d = ', num2str(d)])
            savefig(fighandle4, ...
                "result_"+num2str(d)+"dist_"+num2str(p)+"depth_", ...
                path = "Eigenvalues\20febtest\Result", overwrite = 1)


            max_length = max(cellfun('length',{mv_getCapon.StructTmp.R_eigvals}));
            Reig_axial = zeros(N, max_length);
            Reig_depth = zeros(E, max_length);
            
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
            % Normalize the eigenvalues
            max_eigval = max(max(Reig_axial));
            Reig_axial(:,:) = Reig_axial(:,:)./max_eigval;
            % Replaxe nan values
            Reig_axial(isnan(Reig_axial)) = -0.1;

            max_eigval = max(max(Reig_depth));
            Reig_depth(:,:) = Reig_depth(:,:)./max_eigval;
            % Replaxe nan values
            Reig_depth(isnan(Reig_depth)) = -0.1;
%%
            fighandle1 = figure; 
            fighandle1.Position = [100 100 1000 400];

            subplot(121)
            imagesc(depth_axis*1e3, 1:max_length, db(Reig_axial).');
            ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
            xlabel("Depth [mm]")
            title("Through " + num2str(az) + " degrees")
            colbar = colorbar();
            colbar.Title.String = "Eigenvalue [dB]";
            clim([-60 0])
            xlim([20 40])
            subplot(122)
            imagesc(-20:20, 1:max_length, db(Reig_depth).')
            colbar = colorbar();
            colbar.Title.String = "Eigenvalue [dB]";
            clim([-60 0])
            title("Through 30 mm")
            ylabel("\ $\left| \lambda_i \right|$",'Interpreter','latex')
            xlabel("Width [mm] EGENTLIG IKKE SÅ FIKS SENERE")
%             xlim([0 200])
            
            sgtitle("Eigenvalues sorted, Subarray length: " + num2str(Lelm) + ", d="+num2str(d))

            savefig(fighandle1, ...
                "eigvals_"+num2str(d)+"dist_"+num2str(p)+"depth_", ...
                path = "Eigenvalues\20febtest\Eigenvalues", overwrite = 1)

            %% Gather resolution
%             tmpstruct = uff.beamformed_data(b_data_mv_getCapon);
%             tmpstruct.data = b_data_mv_getCapon.data;


            res_inds = depth_inds:length(depth_axis):length(b_data_mv_getCapon.data);
%             tmpstruct.data(res_inds) = 10;

%             tmpstruct.plot(figure(1234))

            fig_resolution = figure(101);
            fig_resolution.Position = [100 100 1000 400];
            plot(azimuth_axis, db(abs(b_data_mv_getCapon.data(res_inds)./max(b_data_mv_getCapon.data(res_inds)))), 'DisplayName', ['d=',num2str(d)])
            hold on
        end
    end 
end
%%
plot(azimuth_axis, -6*ones(length(azimuth_axis),1), '-r', 'DisplayName','-6dB')
hold off
title("Resolution through 30mm")
xlabel("Azimuth")
ylabel("Depth [mm]")
ylim([-60 5])
legend('Location', 'southeast')
savefig(fig_resolution, ...
    "resolution_"+num2str(d)+"dist_", ...
    path = "Eigenvalues\20febtest\Resolution", overwrite = 1)
            %% Plot condition number
%             CN = zeros(N,E);
%             CN_DL = zeros(N,E);
%             Tse = zeros(N,E);
%             Reigval = zeros(N,E);
%             Reigval_DL = zeros(N,E);
%             for i = 1:N
%                 for j = 1:E
%                     R_tmp = mv_getCapon.StructTmp(i,j).R;
%                     R_DL_tmp = mv_getCapon.StructTmp(i,j).R_DL;
%                     w_tmp = mv_getCapon.StructTmp(i,j).w;
%                     M_new_tmp = mv_getCapon.StructTmp(i,j).M;
%                     [CN(i,j), CN_DL(i,j), Tse(i,j), Reigval(i,j), Reigval_DL(i,j), R_eigval] = parameter_calculations(R_tmp, R_DL_tmp, w_tmp, M_new_tmp);
%                 end
%             end
%             b_data_CondNr                =   uff.beamformed_data(b_data_mv_getCapon);
%             b_data_CondNr_DL             =   uff.beamformed_data(b_data_mv_getCapon);
%             b_data_CondNr.data(:,:)      =   CN(:);
%             b_data_CondNr_DL.data(:,:)   =   CN_DL(:);
%             fighandle2 = figure;
%             b_data_CondNr.plot(fighandle2, ['Condition Number, pos = [', num2str(d), ', +-',num2str(p), ']mm, Lelm= ', num2str(Lelm)], [], 'none')
%             fighandle3 = figure;
%             b_data_CondNr_DL.plot(fighandle3, ['Condition Number, pos = [', num2str(d), ', +-',num2str(p), ']mm, Lelm= ', num2str(Lelm)], [], 'none')
     
%             savefig(fighandle2, ...
%                 "CN_"+num2str(d)+"dist_"+num2str(p)+"depth_"+num2str(Lelm)+"L", ...
%                 path = "Condition Number\First test", overwrite = 1)
%             savefig(fighandle3, ...
%                 "CNDL_"+num2str(d)+"dist_"+num2str(p)+"depth_"+num2str(Lelm)+"L", ...
%                 path = "Condition Number\First test", overwrite = 1)

            %%
%             
%             depth_inds = find(abs(depth_axis-30e-3)<1e-4);
%             depth_inds = depth_inds(round(length(depth_inds)/2));
%             
%             fig_condnr_depth = figure(100);
%             plot(rad2deg(azimuth_axis).', CN(depth_inds,:), 'DisplayName',['L=',num2str(Lelm)])
%             xlabel("Azimuth [degrees]")
%             ylabel("Condition number")
%             hold on
%             fig_condnr_az = figure(101);
%             plot(depth_axis*1e3, CN(:, az_inds), 'DisplayName',['L=',num2str(Lelm)])
%             xlabel("Depth [mm]")
%             ylabel("Condition number")
%             hold on
%         end
%     end
% end

%% 
% az = deg2rad(0); % deg
% az_inds = find(abs(azimuth_axis-az)<0.01);
% az_inds = az_inds(round(length(az_inds)/2));
% 
% depth_inds = find(abs(depth_axis-30e-3)<1e-4);
% depth_inds = depth_inds(round(length(depth_inds)/2));
% 
% fig_condnr_depth = figure(1000);
% plot(rad2deg(azimuth_axis), CN(depth_inds,:), 'DisplayName','L='+num2str(Lelm))
% xlabel("Azimuth [degrees]")
% ylabel("Condition number")
% fig_condnr_az = figure(1001);
% plot(depth_axis*1e3, CN(:, az_inds), 'DisplayName','L='+num2str(Lelm))
% xlabel("Depth [mm]")
% ylabel("Condition number")
% 
% fighandle4 = figure();
% b_data_mv_getCapon.plot(fighandle4, ['getCapon Output, Lelm = ', num2str(Lelm)])