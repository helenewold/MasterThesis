%% analysis_ddist.m
% Script to generate points at a distance d from each other.
clear;
close all;
% if 
%     ispc
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\downloads\Field_II_ver_3_30_windows.tar'));
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\pc\Dokumenter\USTB'));
%     addpath(genpath('\\hume.uio.no\student-u55\helenewo\MasterThesis'));
%     path_ = ['\Figures\Server\',date];
%     path_fig = ['\FiguresFigFormat\Server\',date];
%     basepath = '\\hume.uio.no\student-u55\helenewo\MasterThesis\MasterThesis';
if isunix
    addpath(genpath('/hom/dsb/field'));
    addpath(genpath('/uio/hume/student-u55/helenewo/pc/Dokumenter/USTB'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis'));
    addpath(genpath('/uio/hume/student-u55/helenewo/MasterThesis/datasets'));
    path_ = ['./Figures/',date];
    path_fig = ['./FiguresFigFormat/',date];
    basepath = '/uio/hume/student-u55/helenewo/MasterThesis/MasterThesis';
end
basepath = 'C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis';


p = 30;%[20, 30, 40];
dist = [0, 0.5, 1, 1.5, 2, 5];

degs_separated = 360 * (dist / (2*pi*p));

Lfrac = 1/3;

MLAs = 3;


for MLA = MLAs
    DistResMVAmp = zeros(length(dist), 128*MLA);
    DistResMVPow = zeros(length(dist), 128*MLA);
    DistResDAS = zeros(length(dist), 128*MLA);
    maxDasVals = zeros(length(dist), 1);
    ii=1;

    for degs_i = 1:length(degs_separated)
        degs = degs_separated(degs_i);
        d = dist(degs_i);

%         if d == 0
%             phantom_positions = [0, 0, p*1e-3];
%         else
        phantom_positions(1,:) = [ p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
        phantom_positions(2,:) = [-p*1e-3*sind( degs*0.5 ), 0, p*1e-3*cosd( degs*0.5 )];
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
            
            depth_inds = find(abs(depth_axis-p*1e-3)<1e-4);
            depth_inds = depth_inds(round(length(depth_inds)/2));


%             axises = zeros(N*E,1);
%             axises(depth_inds:N:end) = 1;           % Depth
%             axises(N*az_inds:1:N*(az_inds+1)) = 1;  % Azimuth
%             inds = find(axises == 1);
%             b_data_mv_getCapon.data(inds) = 10;
% 
% 
            fighandle4 = figure();
            b_data_mv_getCapon.plot(fighandle4, ['getCapon Output, d = ', num2str(d)])
%             save_fig(fighandle4, ...
%                 "result_"+num2str(d)+"dist_"+num2str(p)+"depth_"+num2str(MLA)+"MLA", ...
%                 path = "Eigenvalues\"+date+"\Result")


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
            % Normalize the eigenvalues
%             max_eigval = max(max(Reig_axial));
%             Reig_axial(:,:) = Reig_axial(:,:)./max_eigval;
            % Replaxe nan values
            Reig_axial(isnan(Reig_axial)) = -0.1;

%             max_eigval = max(max(Reig_depth));
%             Reig_depth(:,:) = Reig_depth(:,:)./max_eigval;
            % Replaxe nan values
            Reig_depth(isnan(Reig_depth)) = -0.1;
%%
            fighandle1 = figure; 
%             fighandle1.Position = [100 100 750 450];
%             subplot(121)
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

%             save_fig(fighandle1, ['EV_',num2str(d),'dist_azim'], fig_path = fullfile(date,'Eigenvalues', 'Azim'), path = fullfile(date,'Eigenvalues', 'Azim'))
%             save_fig(fighandle2, ['EV_',num2str(d),'dist_depth'], fig_path = fullfile(date,'Eigenvalues', 'Depth'), path = fullfile(date,'Eigenvalues', 'Depth'))
            

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
    end 
%%
maxDas = maxDasVals(end);

% DistRes(ii,:) = DAStmp(depth_inds:N:end);

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
    plot(x_axis_dmm, db(abs(DistResMVAmp(iii,:))./maxDas), 'LineWidth',1.25, 'LineStyle', lt(iii),'DisplayName', string(L(iii)))
    hold on
end
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm] - Amplitude Capon")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([-60 6])
xlim([-6 6])
yticks(-60:6:6)
grid on

%%
figure_res_MVPow = figure(1235);
drawnow;
figure_res_MVPow.Position = [200 100 1000 400];
for iii = 1:length(dist)
    plot(x_axis_dmm, db(abs(DistResMVPow(iii,:))), 'LineWidth',1.25, 'LineStyle', lt(iii),'DisplayName', string(L(iii)))
    hold on
end
hold off

legend('Location', 'northeast')
title("Lateral line through 30[mm] - Power Capon")
xlabel("x [mm]")
ylabel("Amplitude [dB]")
ylim([max(max(db(abs(DistResMVPow))))-60 max(max(db(abs(DistResMVPow))))+6])
xlim([-6 6])
yticks(-84:6:6)
grid on

%%

figure_res_DAS = figure(1236);
drawnow;
figure_res_DAS.Position = [200 100 1000 400];
for iii = 1:length(dist)
    plot(x_axis_dmm, db(abs(DistResDAS(iii,:)./maxDas)), 'LineWidth',1.25, 'LineStyle', lt(iii),'DisplayName', string(L(iii)))
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
path_ = fullfile(date,'Ddist');
% save_fig(figure_res_MVAmp, 'resolution_AmpCapon', path = path_, fig_path =path_)
% save_fig(figure_res_MVPow, 'resolution_PowCapon', path = path_, fig_path =path_)
% save_fig(figure_res_DAS,   'resolution_DAS', path = path_, fig_path =path_)
if ~exist(fullfile(basepath, 'MT SourceCode', 'Data', 'Fullrun', 'SpeckleLevel0', 'ddist')); mkdir(fullfile(basepath, 'MT SourceCode', 'Data', 'Fullrun', 'SpeckleLevel0', 'ddist')); end

save(fullfile(basepath, 'MT SourceCode', 'Data', 'Fullrun', 'SpeckleLevel0', 'ddist', 'ddist_nospeckle.mat'), 'dist', "DistResDAS", 'DistResMVAmp', 'DistResMVPow', 'x_axis_dmm', "maxDas")



%%

dist = [dist(1:5) dist(end)];
nplots = length(dist);

fig_azimuth = figure(2345);
fig_azimuth.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname_az      = "..\Figures\29-Mar-2023\Eigenvalues\Azim\EV_"  +num2str(d)  +   "dist_azim.fig";
    fig=openfig(genname_az);
    pause(0.5)

%     clim([-60 0])
    xlim([25 34])
    ylim([1 10])
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d), ' mm'];
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    if i~=1&&i~=4
        set(ax,'YTickLabel',[]);
    end    
    if i<=3
        set(ax,'XTickLabel',[]);
    end
    
    ax.Parent = tcl;
    ax.Layout.Tile = i;
    set(ax, 'CLim', [-60 5])

    i = i+1;
    close(fig)

end
fig_azimuth;
sgtitle("Eigenvalues of R through axial line at 0 degrees")

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Label.String = "Eigenvalue [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3;
% cb.set('Limits', [-60 0])

tcl.XLabel.String = 'Depth [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 14;
%%


fig_depth = figure(234);
fig_depth.Position = [100 200 800 400];
tcl=tiledlayout(2, 3, 'TileSpacing','tight','Padding','compact');


i = 1;
for d = dist
    genname_az      = "..\Figures\29-Mar-2023\Eigenvalues\Depth\EV_"  +num2str(d)  +   "dist_depth.fig";
    fig=openfig(genname_az);
    pause(0.5)

%     clim([-60 0])
%     xlim([25 34])
    ylim([1 10])
    colorbar off
    ax = gca;
    ax.Title.String = ['d= ', num2str(d), ' mm'];
    ax.XLabel.Visible = 'off';
    ax.YLabel.Visible = 'off';
    if i~=1&&i~=4
        set(ax,'YTickLabel',[]);
    end    
    if i<=3
        set(ax,'XTickLabel',[]);
    end
    
    ax.Parent = tcl;
    ax.Layout.Tile = i;
    set(ax, 'CLim', [-60 5])

    i = i+1;
    close(fig)

end
fig_depth;
sgtitle("Eigenvalues of R at lateral line through 30[mm]")

cb = colorbar();
cb.Layout.Tile = 'east';
cb.Label.String = "Eigenvalue [dB]";
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3;

tcl.XLabel.String = 'x [mm]';
tcl.XLabel.FontSize = 16;

tcl.YLabel.Interpreter = 'latex';
tcl.YLabel.String = "$\left| \lambda_i \right|$";
tcl.YLabel.FontSize = 14;

%% Save
basepath = "";
path_ = fullfile(date,'Eigenvalues');
% save_fig(fig_depth, 'EV_ddist_depth', path = path_, fig_path = path_)%, base_path = basepath)
% save_fig(fig_azimuth, 'EV_ddist_azimuth', path = path_, fig_path = path_)%, base_path = basepath)
% path =  fullfile(basepath, 'Figures', date,'Ddist', 'Eigenvalues');
% save_fig(fig_depth, 'EV_Lelm_depth_MLA5', path = path, fig_path = path, base_path = '')%, base_path = basepath)
% save_fig(fig_azimuth, 'EV_Lelm_azimuth_MLA5', path = path, fig_path = path, base_path = '')%, base_path = basepath)


end


%%
% plot(x_axis_dmm, -6*ones(length(x_axis_dmm),1), '-r', 'DisplayName','-6dB')
% hold off
% title("Resolution through 30mm")
% xlabel("Width [mm]")
% ylabel("Amplitude [dB]")
% ylim([-60 5])
% xlim([x_axis_dmm(1) x_axis_dmm(end)])
% legend('Location', 'southeast')
% savefig(fig_resolution, ...
%     "resolution_", ...
%     path = "Eigenvalues\"+date+"\Resolution")
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