function K_samples = TimeAverageCalculation(scan, channel_data, K_in_lambda)

    if isa(scan,'uff.linear_scan')
        z_in_lambda = h.scan(1).z_axis./h.channel_data.lambda;
        z_in_lambda = z_in_lambda - z_in_lambda(1);
        [~,samples] = min(abs(z_in_lambda-K_in_lambda));
    elseif isa(scan,'uff.sector_scan')
        z_in_lambda = scan(1).depth_axis./channel_data.lambda;
        z_in_lambda = z_in_lambda - z_in_lambda(1);
        [~,samples] = min(abs(z_in_lambda-K_in_lambda));
    end
    
    if mod(round(samples),2)    % Check if odd
        K_samples = round(samples);
    else
        K_samples = round(samples)+1;
    end
end