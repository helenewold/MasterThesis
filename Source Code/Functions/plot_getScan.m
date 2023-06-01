function plot_getScan(b_data, azimuth, depth, options)
    arguments
        b_data
        azimuth
        depth
        options.title = ""
        options.xlabel = "Width [mm]"
        options.ylabel = "Depth [mm]"
        options.dB = 1
        options.dyn_range = 50
    end

    [getScan, Xs, Zs] = getScanConvertedImage(b_data.data, azimuth, ...
                                                        depth, length(azimuth), ...
                                                        length(depth));%, sizeX, sizeZ, interpolationMethod);
    inds = isnan(getScan);

    if options.dB == 1
        getScan(inds) = 1;
        imagesc(Xs*1e3, Zs*1e3, db(abs(getScan)))
    elseif options.dB == 2
        getScan(inds) = 0;
        imagesc(Xs*1e3, Zs*1e3, log10(getScan))
    else
        getScan(inds) = 0;
        imagesc(Xs*1e3, Zs*1e3, getScan)
    end
    
    dyn_range = [max(max(db(abs(getScan(:,:)))))-options.dyn_range, max(max(db(abs(getScan(:,:)))))];
%     dyn_range = [-options.dyn_range, 0];

    xlabel(options.xlabel)
    ylabel(options.ylabel)
    title(options.title)
    colormap(gray)
    hColourbar = colorbar();
    if options.dB == 1
        hColourbar.Label.String = "dB";
        hColourbar.Label.Position(1) = 2.75;
        hColourbar.Label.Rotation = 270;
        caxis(dyn_range)
    end
end