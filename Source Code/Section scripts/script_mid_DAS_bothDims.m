azimuth_MLA = linspace(azimuth_axis(1), azimuth_axis(end), length(azimuth_axis)*MLA)';
scan_MLA = uff.sector_scan('azimuth_axis',linspace(azimuth_axis(1), azimuth_axis(end), length(azimuth_axis)*MLA)' ,'depth_axis',depth_axis);

mid_DAS = midprocess.das();
mid_DAS.channel_data = channel_data;
mid_DAS.dimension = dimension.both();

mid_DAS.receive_apodization.window = receive_window;

if exist('f_number', 'var')
    mid_DAS.receive_apodization.f_number = f_number;
end

mid_DAS.pw_margin = pw_margin;
mid_DAS.spherical_transmit_delay_model = spherical_transmit_delay_model_;
mid_DAS.transmit_apodization.window = transmit_window;
mid_DAS.transmit_apodization.MLA = MLA;
mid_DAS.transmit_apodization.MLA_overlap = MLA_overlap;
mid_DAS.scan = scan_MLA;

b_data_DAS = mid_DAS.go();