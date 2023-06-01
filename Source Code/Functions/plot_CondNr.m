function plot_CondNr(b_data_CondNr, azimuth, depth, options)
    arguments
        b_data_CondNr
        azimuth
        depth
        options.title = ""
        options.xlabel = "Width [mm]"
        options.ylabel = "Depth [mm]"
        options.ax
        options.dB = 'none' %'none', 'log10', 'section', '%'
    end

    tmp = b_data_CondNr.data;
    if strcmp(options.dB, 'section')
        tmp(tmp < 1e2 ) = 2;
        tmp(b_data_CondNr.data>1e2 & b_data_CondNr.data<1e7) = 4;
        tmp(b_data_CondNr.data>1e7 & b_data_CondNr.data<1e16) = 6;
        tmp(b_data_CondNr.data > 1e16) = 8;
    end

    
    [getScan, Xs, Zs] = getScanConvertedImage(tmp, azimuth, ...
                                                depth, length(azimuth), ...
                                                length(depth));
    getScan(isnan(getScan)) = 0;

    maxval = max(tmp);

    custom_color(1,:)   =   repmat([1 1 1],                 1, 1);     % hvit bakgrunn
    custom_color(2:3,:) =   repmat([0 1 0],                 2, 1);     % Grønn
    if maxval == 4
        custom_color(4,:) =   repmat([1 1 0],                 1, 1);     % gul
    elseif maxval == 6
        custom_color(4:5,:) =   repmat([1 1 0],                 2, 1);     % gul
        custom_color(6,:) =   repmat([0.929 0.694 0.125],     1, 1);     % oransje
    elseif maxval == 8
        custom_color(4:5,:) =   repmat([1 1 0],                 2, 1);     % gul
        custom_color(6:7,:) =   repmat([0.929 0.694 0.125],     2, 1);     % oransje
        custom_color(8,:)   =   repmat([1 0 0],                 1, 1);     % rød
    end
    cMap = interp1([0;1],[0 1 0; 1 0 0],linspace(0,1,256));


    if strcmp(options.dB, 'none')
        imagesc(Xs*1e3, Zs*1e3, log10(abs(getScan)))
        if isfield(options, 'ax')
            colormap(options.ax, cMap)
        else
            colormap(cMap)
        end
        colbar              =   colorbar();
    elseif strcmp(options.dB, 'section')
        imagesc(Xs*1e3, Zs*1e3, getScan)

        if isfield(options, 'ax')
            colormap(options.ax, custom_color)
        else
            colormap(custom_color)
        end
        colbar              =   colorbar();
        colbar.Limits       =   [1 maxval+1];
        colbar.Ticks        =   0:2:(maxval+1);
        colbar.TickLabels   =   ["Background", "Stable", " Semi stable", "Unstable", "Uninvertable"];

    elseif strcmp(options.dB, 'log10')
        imagesc(Xs*1e3, Zs*1e3, log10(abs(getScan)))
        if isfield(options, 'ax')
            colormap(options.ax, cMap)
        else
            colormap(cMap)
        end
        colbar              =   colorbar();
    elseif strcmp(options.dB, '%')
        % Prosent av egenverdier
        imagesc(Xs*1e3, Zs*1e3, getScan)
        if isfield(options, 'ax')
            colormap(options.ax, cMap)
        else
            colormap(cMap)
        end
        colbar              =   colorbar();
    end

    xlabel(options.xlabel)
    ylabel(options.ylabel)
    title(options.title)


end