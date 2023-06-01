mv_getCapon = postprocess.capon_MV_getCapon();
mv_getCapon.dimension = dimension.receive;

mv_getCapon.transmit_apodization= mid.transmit_apodization;
mv_getCapon.receive_apodization = mid.receive_apodization;
mv_getCapon.scan = scan_MLA;

mv_getCapon.channel_data = channel_data;

mv_getCapon.K_in_lambda = K_in_lambda;
mv_getCapon.L_elements = Lelm;


if exist("Lelm_set", 'var')
    mv_getCapon.L_elements_set = Lelm_set;
else
    mv_getCapon.L_elements_set = 0;
end


if exist("USTB_normalization", 'var')
    mv_getCapon.USTB_normalization = USTB_normalization;
else
    mv_getCapon.USTB_normalization = 1;
end

if exist("calc_param", 'var')
    mv_getCapon.calc_param = calc_param;
% else
%     mv_getCapon.calc_param = 0;
end

% mv_getCapon.fig_handle = fighandle;


mv_getCapon.regCoef = regCoeff;

mv_getCapon.input = b_data_MLA_delayed;

[rx_apodization, tx_apodization] = create_apod_matrix(mv_getCapon);

% b_data_mv_getCapon = mv_getCapon.go();
