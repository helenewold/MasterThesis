mv_MLA = postprocess.capon_minimum_variance();
mv_MLA.dimension = dimension.receive;

mv_MLA.transmit_apodization= mid.transmit_apodization;
mv_MLA.receive_apodization = mid.receive_apodization;
mv_MLA.scan = scan_MLA;

mv_MLA.channel_data = channel_data;

mv_MLA.K_in_lambda = K_in_lambda;
mv_MLA.L_elements = Lelm;

mv_MLA.regCoef = regCoeff;



mv_MLA.input = b_data_MLA_delayed;
b_data_mv_MLA = mv_MLA.go();
