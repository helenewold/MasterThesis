function [rx_apodization, tx_apodization] = create_apod_matrix(PrepStruct)

    % Initialize
    % check if we have information about apodization
    rx_apodization=ones([PrepStruct.input(1).N_pixels,PrepStruct.input.N_channels]);
    tx_apodization=ones([PrepStruct.input(1).N_pixels,PrepStruct.input.N_waves]);
    %%setter rx_apodization
    if ~isempty(PrepStruct.transmit_apodization)&~isempty(PrepStruct.receive_apodization)&~isempty(PrepStruct.channel_data.probe)
        % receive
        if PrepStruct.input.N_channels > 1
            PrepStruct.receive_apodization.probe=PrepStruct.channel_data.probe;
            rx_apodization=PrepStruct.receive_apodization.data();
        end
        
        % transmit
        if PrepStruct.input.N_waves > 1
            PrepStruct.transmit_apodization.sequence = PrepStruct.channel_data.sequence;
            PrepStruct.transmit_apodization.probe=PrepStruct.channel_data.probe;
            tx_apodization=PrepStruct.transmit_apodization.data();
        end
    else
        warning('Missing probe and apodization data; full aperture is assumed.');
    end
end