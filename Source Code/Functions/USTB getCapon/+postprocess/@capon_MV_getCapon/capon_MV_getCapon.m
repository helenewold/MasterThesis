classdef capon_MV_getCapon < postprocess
    %   CAPON_MINIMUM_VARIANCE getCapon
    %
    %             Beamform a ultrasound image using the Capon beamformer
    %
    %   implementers: Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %                 Helene Wold
    %
    %   $Last updated: 2022$
    
    %% constructor
    methods (Access = public)
        function h=capon_MV_getCapon()
            h.name='Minimum Variance MATLAB';
            h.reference= 'Benefits';
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>'};
            h.version='v1.0.0';
        end
    end
    
    %% Additional properties
    properties
        active_element_criterium=0.16;                % value to decide whether an element is used or not
        L_elements                                    % subarray size
        
        K_in_lambda                                   % temporal averaging factor
        regCoef                                       % regularization factor
        doForwardBackward = 0;                        % forward backward averaging
        dimension = dimension.both;                   % dimension class that specifies whether the process will run only on transmit, receive, or both.
        channel_data                                  % Channel data 
        scan
        L_elements_set = 0;
        indsI=0
        indsJ = 0
        V = []
        usePinv = false

        USTB_normalization = 1;
        amp = 1;
        
        calc_param = 0;
        K_samples
        
        imPower
        CN
        CN_DL
        Tse
        Reigval
        Reigval_DL
        CosineSimilarity
        
        subarrsize
        StructTmp
%         fig_handle = figure();
    end
    
    methods
        function [output, imPower]=go(h)       
            % check if we can skip calculation
            if h.check_hash()
                output = h.output; 
                imPower = h.imPower;
                return;
            end  
            
            % check dimensions
            if (h.dimension==dimension.receive) && (h.input.N_channels<2)
                error('Not enough channels to compute factor');
            end
            if (h.dimension==dimension.transmit) && (h.input.N_waves<2)
                error('Not enough waves to compute factor');
            end
            if (h.dimension==dimension.both) 
                if (h.input.N_channels<2)&&(h.input.N_waves>1)
                    warning('Not enough channels to compute factor. Changing dimension to dimension.transmit');
                    h.dimension = dimension.transmit;
                elseif (h.input.N_waves<2)&&(h.input.N_channels>1)
                    warning('Not enough waves to compute factor. Changing dimension to dimension.receive');
                    h.dimension = dimension.receive;
                elseif (h.input.N_waves<2)&&(h.input.N_channels<2)
                    error('Not enough waves and channels to compute factor');
                end
            end
            
            % check if we have information about apodization
            rx_apodization=ones([h.input(1).N_pixels,h.input.N_channels]);
            tx_apodization=ones([h.input(1).N_pixels,h.input.N_waves]);
            if ~isempty(h.transmit_apodization)&~isempty(h.receive_apodization)&~isempty(h.channel_data.probe)
                % receive
                if h.input.N_channels > 1
                    h.receive_apodization.probe=h.channel_data.probe;
                    rx_apodization=h.receive_apodization.data();
                end
                
                % transmit
                if h.input.N_waves > 1
                    h.transmit_apodization.sequence = h.channel_data.sequence;
                    h.transmit_apodization.probe=h.channel_data.probe;
                    tx_apodization=h.transmit_apodization.data();
                end
            else
                warning('Missing probe and apodization data; full aperture is assumed.');
            end  

            % declare output structure
            h.output=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data
            
            switch h.dimension
                case dimension.both
                    str = ['You are trying to run the Capon minimum variance beamformer on both dimensions simultaneously. ',...
                           'This is to my knowledge not been done in the litterature before, and might not make sense. ',...
                           'I also takes forever...'];
                    warning(str);

                    % auxiliary data
                    aux_data=zeros(h.input.N_pixels,1,1,h.input.N_frames);
                    aux_data_Power=zeros(h.input.N_pixels,1,1,h.input.N_frames);
                    for n_frame = 1:h.input.N_frames
                        apod_matrix = zeros(size(tx_apodization,1),h.input.N_waves*h.input.N_channels);
                        for i = 1:h.input.N_waves
                            apod_matrix(:,1+(i-1)*h.input.N_channels:h.input.N_channels*i) = tx_apodization(:,i).*rx_apodization;
                        end
                        apod_matrix = reshape(apod_matrix,h.input(1).scan.N_z_axis,h.input(1).scan.N_x_axis,h.input.N_waves*h.input.N_channels);
%                     
                        % Apodization matrix indicating active elements
                        %apod_matrix = reshape(bsxfun(@times,tx_apodization,rx_apodization(n_channel)),h.input(1).scan.N_z_axis,h.input(1).scan.N_x_axis,h.input.N_waves);
                        data_cube = reshape(h.input.data(:,:,:,n_frame),h.input(1).scan.N_z_axis,h.input(1).scan.N_x_axis,h.input.N_channels*h.input.N_waves);
                        [imAmplitude, imPower] = get_Capon(data_cube, apod_matrix, h, ['1/1']);

                        aux_data(:,1,1,n_frame) = imAmplitude(:); 
                        aux_data_Power(:,1,1,n_frame) = imPower(:); 
                    end
                    h.output.data = aux_data;
                case dimension.transmit
                    % auxiliary data
                    aux_data=zeros(h.input.N_pixels,h.input.N_channels,1,h.input.N_frames);
                    for n_frame = 1:h.input.N_frames
                        for n_channel = 1:h.input.N_channels
                            % Apodization matrix indicating active elements
                            apod_matrix = reshape(bsxfun(@times,tx_apodization,rx_apodization(n_channel)),h.input(1).scan.N_z_axis,h.input(1).scan.N_x_axis,h.input.N_waves);
                            data_cube = reshape(h.input.data(:,n_channel,:,n_frame),h.input(1).scan.N_z_axis,h.input(1).scan.N_x_axis,h.input.N_waves);
                            [imAmplitude, imPower] = get_Capon(data_cube, apod_matrix, h, [num2str(n_channel),'/',num2str(h.input.N_channels)]);

                            aux_data(:,n_channel,:,n_frame) = imAmplitude(:);
                            aux_data_Power(:,1,1,n_frame) = imPower(:); 
                        end
                    end
                    h.output.data = aux_data;
                case dimension.receive
                    % auxiliary data
                    aux_data=zeros(h.input.N_pixels,1,h.input.N_waves,h.input.N_frames);
                    for n_frame = 1:h.input.N_frames
                        for n_wave = 1:h.input.N_waves
                            % Apodization matrix indicating active elements
                            if isa(h.scan,'uff.linear_scan')
                                apod_matrix = reshape(bsxfun(@times,tx_apodization(:,n_wave),rx_apodization),h.input(1).scan.N_z_axis,h.input(1).scan.N_x_axis,h.input.N_channels);
                                data_cube = reshape(h.input.data(:,:,n_wave,n_frame),h.input(1).scan.N_z_axis,h.input(1).scan.N_x_axis,h.input.N_channels);
                                
                            elseif isa(h.scan,'uff.sector_scan')
                                apod_matrix = reshape(bsxfun(@times,tx_apodization(:,n_wave),rx_apodization),h.input(1).scan.N_depth_axis,h.input(1).scan.N_azimuth_axis,h.input.N_channels);
                                data_cube = reshape(h.input.data(:,:,n_wave,n_frame),h.input(1).scan.N_depth_axis,h.input(1).scan.N_azimuth_axis,h.input.N_channels);
                                
                            end
                                                       %A hack to set non active elements to zero for the
                            %alpinion scanner FI who only use 64 active
                            %elements
                            if ~isempty(h.channel_data.N_active_elements) && sum(h.channel_data.N_active_elements ~= h.channel_data.N_elements)
                                apod_matrix(abs(data_cube)<eps) = 0;
                            end
                            [imAmplitude, imPower] = get_Capon(h, data_cube, apod_matrix, [num2str(n_wave),'/',num2str(h.input.N_waves)]);
                            aux_data(:,1,n_wave,n_frame) = imAmplitude(:);
                            aux_data_Power(:,1,1,n_frame) = imPower(:); 
                        end
                    end
                    inds = find(isnan(aux_data));
                    aux_data(inds) = 0;
                    inds = find(isnan(aux_data_Power));
                    aux_data_Power(inds) = 0;

                    h.imPower = aux_data_Power;
                    h.output.data = aux_data;
                otherwise
                    error('Unknown dimension mode; check HELP dimension');
            end
            
            % pass reference
            output = h.output;
            imPower = h.imPower;
            
            % update hash
            h.save_hash();
        end
        
        function h=set.K_in_lambda(h,K_in_lambda)
            assert(~isempty(h.scan),'You need to set the scan.')
            assert(~isempty(h.channel_data),'You need to set the channel_data.')
          
            if isa(h.scan,'uff.linear_scan')
                h.K_in_lambda = K_in_lambda;
                z_in_lambda = h.scan(1).z_axis./h.channel_data.lambda;
                z_in_lambda = z_in_lambda - z_in_lambda(1);
                [~,samples] = min(abs(z_in_lambda-h.K_in_lambda));
            elseif isa(h.scan,'uff.sector_scan')
                h.K_in_lambda = K_in_lambda;
                z_in_lambda = h.scan(1).depth_axis./h.channel_data.lambda;
                z_in_lambda = z_in_lambda - z_in_lambda(1);
                [~,samples] = min(abs(z_in_lambda-h.K_in_lambda));
            end
            
            if mod(round(samples),2)    % Check if odd
                h.K_samples = round(samples);
            else
                h.K_samples = round(samples)+1;
            end
        end

    end
end



