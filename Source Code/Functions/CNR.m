function CNR_value = CNR(pixels_ROI, pixels_B)
    % takes in pixel in amplitude
    pixels_ROI = abs(pixels_ROI(:));
    pixels_B = abs(pixels_B(:));
    
    u_ROI = mean( pixels_ROI.^2 );
    u_B   = mean( pixels_B.^2 );
    var_ROI = var( pixels_ROI.^2 );
    var_B   = var( pixels_B.^2 );
    

    CNR_value = abs(u_ROI - u_B)  / sqrt(var_ROI + var_B + eps); 
end

