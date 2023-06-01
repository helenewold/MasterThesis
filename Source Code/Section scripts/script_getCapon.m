tmp_data=zeros(post_getCapon.input.N_pixels,1,post_getCapon.input.N_waves,post_getCapon.input.N_frames);


for n_frame = 1:post_getCapon.input.N_frames
    for n_wave = 1:post_getCapon.input.N_waves
        % Apodization matrix indicating active elements
        if isa(post_getCapon.scan,'uff.linear_scan')
            apod_matrix = reshape(bsxfun(@times,tx_apodization(:,n_wave),rx_apodization), ...
                post_getCapon.input(1).scan.N_z_axis, ...
                post_getCapon.input(1).scan.N_x_axis, ...
                post_getCapon.input.N_channels);
            data_cube = reshape(post_getCapon.input.data(:,:,n_wave,n_frame), ...
                post_getCapon.input(1).scan.N_z_axis, ...
                post_getCapon.input(1).scan.N_x_axis, ...
                post_getCapon.input.N_channels);
        elseif isa(post_getCapon.scan,'uff.sector_scan')
            apod_matrix = reshape(bsxfun(@times,tx_apodization(:,n_wave),rx_apodization), ...
                post_getCapon.input(1).scan.N_depth_axis, ...
                post_getCapon.input(1).scan.N_azimuth_axis, ...
                post_getCapon.input.N_channels);
            data_cube = reshape(post_getCapon.input.data(:,:,n_wave,n_frame), ...
                post_getCapon.input(1).scan.N_depth_axis, ...
                post_getCapon.input(1).scan.N_azimuth_axis, ...
                post_getCapon.input.N_channels);
        end
        % A hack to set non active elements to zero for the
        % alpinion scanner FI who only use 64 active elements
        if ~isempty(post_getCapon.channel_data.N_active_elements) && sum(post_getCapon.channel_data.N_active_elements ~= post_getCapon.channel_data.N_elements)
            apod_matrix(abs(data_cube)<eps) = 0;
        end

        [image, ~] = getCapon(data_cube, apod_matrix = apod_matrix, verbose = 1, regCoef = 1/100, L = round(channel_data.probe.N/3), nTimeAverage = K_samples);
        tmp_data(:,1,n_wave,n_frame) = image(:);
    end
end

post_getCapon.data = tmp_data;

post_getCapon.N_channels = size(tmp_data,2);
post_getCapon.N_pixels = size(tmp_data,1);
post_getCapon.N_waves = post_getCapon.input.N_waves;
post_getCapon.N_frames = post_getCapon.input.N_frames;