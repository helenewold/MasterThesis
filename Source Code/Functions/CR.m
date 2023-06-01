function CR_dB_value = CR(pixels_ROI, pixels_B)
    % takes in pixels in amplitude
    pixels_ROI = pixels_ROI(:);
    pixels_B = pixels_B(:);
    
    CR_value = mean( pixels_ROI.^2 ) / mean( pixels_B.^2 + eps);
    CR_dB_value = 10*log10( CR_value );
end

