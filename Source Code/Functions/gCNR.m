function [gCNR_value] = gCNR(region1, region2, type, bins, return_plot, label)
    %% Calculate gCNR from two pixel populations (region 1 and 2).
    % region1 = pixels in absolute value 
    % region2 = pixels in absolute value 
    % type = 'kde' (default) or 'histogram'
    % return_plot = true or false

    arguments
        region1
        region2
        type = 'kde';
        bins = 100;
        return_plot = false;
        label = [];
    end
    
    % RENAME to region1 = ROI, region2 = B, like:
    %    pixels_ROI = pixels_ROI(:);
    %    pixels_B = pixels_B(:);
    
    
    %% Check and initialize
    if isempty(type); type = 'kde'; end;
    if isempty(bins); bins = 100; end;
    
    region1 = abs(region1(:));
    region2 = abs(region2(:));
    
    PM = max([region1; region2]);
    region1_dB = db(region1/PM +eps); % problem if it is actually 0 !!!
    region2_dB = db(region2/PM +eps);    
    
    % Bins, limit minimum at -150 dB
    min_xs = max( min( [region1_dB; region2_dB] ), -150);
    max_xs =      max( [region1_dB; region2_dB] );
    xs = linspace( min_xs, max_xs, bins);
    
    %% Compute KDE or histogram
    if strcmp(type, 'kde')
        [pdf_i] = ksdensity(region1_dB, xs);
        [pdf_o] = ksdensity(region2_dB, xs);
        OVL  = sum(min([pdf_i ./sum(pdf_i) ; pdf_o ./sum(pdf_o)] ));
        gCNR_value  = 1 - OVL;
    else strcmp(type, 'histogram') | strcmp(type, 'hist');
        [pdf_i] = hist( region1_dB, xs);
        [pdf_o] = hist( region2_dB, xs);
        OVL = sum(min([pdf_i./sum(pdf_i); pdf_o./sum(pdf_o)]));
        gCNR_value = 1 - OVL;
    end
    
    %% Plot
    if return_plot == true
        figure();
        plot(xs, pdf_i./sum(pdf_i),'r-', 'linewidth',2); hold on; grid on;
        plot(xs, pdf_o./sum(pdf_o),'b-', 'linewidth',2);
        hh=area(xs,min([pdf_i./sum(pdf_i); pdf_o./sum(pdf_o)]), 'LineStyle','none');
        hh.FaceColor = [0.6 0.6 0.6];
        xlabel( 'dB' );
        ylabel('Probability');
        legend('p_i','p_o','OVL');
        title(label)
    end
end



