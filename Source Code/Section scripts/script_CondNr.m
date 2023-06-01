b_data_CondNr_Lfrac                   =   uff.beamformed_data(b_data_mv_getCapon);
b_data_CondNr_after_Lfrac             =   uff.beamformed_data(b_data_mv_getCapon);

b_data_CondNr_Lfrac.data(:,:)         =   mv_getCapon.CN;
b_data_CondNr_after_Lfrac.data(:,:)   =   mv_getCapon.CN_after;


fig_CondNr_Lfrac = figure(150);
fig_CondNr_Lfrac.WindowState = 'maximized';
sgtitle(['Condition numbers, L = M*' , num2str(L_frac)])

subplot(121)
plot_CondNr(b_data_CondNr_Lfrac, azimuth_axis, depth_axis, title = "Before diagonal loading")
subplot(122)
plot_CondNr(b_data_CondNr_after_Lfrac, azimuth_axis, depth_axis, title = "After diagonal loading")

fig_CondNrLog10_Lfrac = figure(200);
fig_CondNrLog10_Lfrac.WindowState = 'maximized';
sgtitle(['Condition numbers, log10, L = M*' , num2str(L_frac)])
subplot(121)
plot_CondNr(b_data_CondNr_Lfrac, azimuth_axis, depth_axis, title = "Before diagonal loading", dB = 2)
subplot(122)
plot_CondNr(b_data_CondNr_after_Lfrac, azimuth_axis, depth_axis, title = "After diagonal loading", dB = 2)
