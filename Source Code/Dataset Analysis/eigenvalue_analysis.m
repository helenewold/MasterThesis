%% generate_2pkt_ddist.m
% Script to generate points at a distance d from each other.
clear;
close all;

path = '..\Field II datasets';
file_name = 'd_dist_3x2points.uff';

% path = 'Datasets';
% file_name = 'P4_FI_121444_45mm_focus.uff';
filename = fullfile(path, file_name);

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

save_fig = 0;

%% Read data
script_ReadData

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
% scan=uff.linear_scan('azimuth_axis',azimuth_axis,'depth_axis',depth_axis);


script_mid_DAS_MLA

%% DAS begge dimensjoner - til bruk for Histogram Matching
% b_data_DAS
script_mid_DAS_bothDims


%% USTB postprocess capon minimum variance
script_post_capon_minvar

%% ustb getcapon postprocess
% mv_getCapon
Lelm_set=0;
calc_param = 0;
% fighandle = figure(100);
script_post_getCapon
b_data_mv_getCapon = mv_getCapon.go();

%% Histogram matching
[img_out, b_data_out] = tools.histogram_match(b_data_DAS,b_data_mv_getCapon);


%% Plotting
fig1 = figure(1);
b_data_DAS.plot(fig1, "DAS output");

fig2 = figure(2);
b_data_mv_getCapon.plot(fig2, "getCapon output");

fig3 = figure(3);
b_data_mv_MLA.plot(fig3, "MV Capon output");

% Histogram matching
fig4 = figure(4);
b_data_out.plot(fig4, "Histogram matching (DAS - getCapon)", 60, 'none');
caxis([-60 0])

%% Eigenvalues plotting
% Starts by closing old figures to avoid plotting points on top of each
% other
if exist('fig6', 'var'); close(fig6); clear(fig6); end
if exist('fig7', 'var'); close(fig7); clear(fig7); end

% Decide which axial degree to look at
az = 0; % Degrees
az_inds = find(abs(azimuth_axis-deg2rad(az))<0.01);
az_inds = az_inds(round(length(az_inds)/2));

fig5 = figure(5);
tot = length(b_data_mv_getCapon.data);
depth_testing = tot*az_inds/length(azimuth_axis):1:tot*(az_inds+1)/length(azimuth_axis);

tst = uff.beamformed_data(b_data_mv_getCapon);
tst.data(depth_testing) = 10;

tst.plot(fig5, "Output visualizing eigenvalue axis");

fig6 = figure(6);
for i = 1:length(depth_axis)
    for j = az_inds
        Rtmp = mv_getCapon.StructTmp(i,j).R;
        eigvals = eig(Rtmp);
        plot(depth_axis(i)*ones(length(eigvals),1)*1e3, eigvals, 'o')
        hold on
    end
end
hold off
ylabel("Eigenvalues")
xlabel("Depth [mm]")
title("Eigenvalues along " + num2str(az) + " degrees")

fig7 = figure(7);
for i = 1:length(depth_axis)
    for j = az_inds
        Rtmp = mv_getCapon.StructTmp(i,j).R;
        eigvals = eig(Rtmp);
        plot(depth_axis(i)*ones(length(eigvals),1)*1e3, db(eigvals), 'o')
        hold on
    end
end
hold off
ylabel("Eigenvalues [dB]")
xlabel("Depth [mm]")
title("Eigenvalues along " + num2str(az) + " degrees in dB")

%% (Delvis testing av ny struktur)
% Krever mye for å få endra alle kodene til denne nye måten, beholder
% gammel metode inntil videre for ikke å ødelegge noe
N = mv_getCapon.scan.N_depth_axis;
E = mv_getCapon.scan.N_azimuth_axis;
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


%% Plotting
b_data_CondNr                =   uff.beamformed_data(b_data_mv_getCapon);
b_data_CondNr_DL             =   uff.beamformed_data(b_data_mv_getCapon);
b_data_CondNr.data(:,:)      =   CN(:);
b_data_CondNr_DL.data(:,:)   =   CN_DL(:);


b_data_Reigval               =   uff.beamformed_data(b_data_mv_getCapon);
b_data_Reigval_DL            =   uff.beamformed_data(b_data_mv_getCapon);
b_data_Reigval.data(:,:)     =   Reigval(:);
b_data_Reigval_DL.data(:,:)  =   Reigval_DL(:);


b_data_CondNr.plot(figure(10)       , 'Condition Number'                      , [], 'none')
b_data_CondNr_DL.plot(figure(11)    , 'Condition Number after DL'             , [], 'none')
b_data_Reigval.plot(figure(12)      , '% of eigenvalues below 1e-10'          , [], 'none')
b_data_Reigval_DL.plot(figure(13)   , '% of eigenvalues below 1e-10 after DL' , [], 'none')
%% Condition number along chosen axis
fig20 = figure(20);
subplot(211)
plot(depth_axis*1e3, CN(:,az_inds))
xlabel('Depth [mm]')
ylabel('Condition number value')
title('Before diagonal loading')
subplot(212)
plot(depth_axis*1e3, CN_DL(:,az_inds))
xlabel('Depth [mm]')
ylabel('Condition number value')
title('After diagonal loading')
%%
max_length = max(max(cellfun('length',{mv_getCapon.StructTmp.R_eigvals})));
Reig_axial = zeros(N, max_length);

for i = 1:length(depth_axis)
    for j = az_inds
        len = cellfun('length',{mv_getCapon.StructTmp(i,j).R_eigvals});
        Reig_axial(i, 1:len) = sort(mv_getCapon.StructTmp(i,j).R_eigvals(:), 'descend');
        Reig_axial(i,(len+1):end) = -0.1;
    end
end


figure(100)
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

sgtitle("Eigenvalues sorted, NaN showed as value -0.1 [-20dB]")





