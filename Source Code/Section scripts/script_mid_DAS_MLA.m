azimuth_MLA = linspace(azimuth_axis(1), azimuth_axis(end), length(azimuth_axis)*MLA)';
scan_MLA = uff.sector_scan('azimuth_axis',linspace(azimuth_axis(1), azimuth_axis(end), length(azimuth_axis)*MLA)' ,'depth_axis',depth_axis);

mid = midprocess.das();
mid.channel_data = channel_data;
mid.dimension = dimension.transmit();

mid.receive_apodization.window = receive_window;

if exist('f_number', 'var')
    mid.receive_apodization.f_number = f_number;
end

mid.pw_margin = pw_margin;
mid.spherical_transmit_delay_model = spherical_transmit_delay_model_;
mid.transmit_apodization.window = transmit_window;
mid.transmit_apodization.MLA = MLA;
mid.transmit_apodization.MLA_overlap = MLA_overlap;
mid.scan = scan_MLA;

b_data_MLA_delayed = mid.go();