function data_cube = USTB_reshape(scan, data, channel_data, options)
    arguments
        scan
        data
        channel_data
        options.n_frame
    end
    if isfield(options, 'n_frame'); n_frame = options.n_frame; else; n_frame = length(data.data(1,1,1,:)); end
    
    % Reshaper for å kunne bruke getCapon
    if isa(scan,'uff.linear_scan')
        data_cube = reshape(data.data(:,:,1,n_frame),scan.N_z_axis,scan.N_x_axis,channel_data.N_channels);
    elseif isa(scan,'uff.sector_scan')
        data_cube = reshape(data.data(:,:,1,n_frame),scan.N_depth_axis,scan.N_azimuth_axis,channel_data.N_channels);
    end
end