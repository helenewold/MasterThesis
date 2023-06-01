function save_fig(fig, filename, options)
    arguments
        fig
        filename
        options.path = date
        options.type = '.pdf'
        options.base_path = fullfile('..','Figures')
        options.fig_path = fullfile('..','Figures', 'Figformat')
        options.overwrite
    end
    if exist('C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Source Code','dir')==7
        cd("C:\Users\Helene\Documents\Elektronikk, informatikk og teknologi\Master Thesis\MT SourceCode\Source Code")
    end
    % Create folder if it does not exist
    if ~exist(fullfile(options.base_path, options.path), 'dir');     mkdir(fullfile(options.base_path, options.path));     end
    if ~exist(fullfile(options.base_path, options.fig_path), 'dir'); mkdir(fullfile(options.base_path, options.fig_path)); end
    
    type = options.type;
    full_path = fullfile(options.base_path, options.path, append(filename, type));

    % Check if file already exists.
    if isfile(full_path)
        if ~isfield(options, 'overwrite') || ~options.overwrite
            fignr = 1;
            while isfile(fullfile(options.base_path, options.path, append(filename ,  num2str(fignr) , type)))
                fignr = fignr+1;
            end
            new_name = fullfile(options.base_path, options.path, append(filename ,  num2str(fignr) , type));
            exportgraphics(fig, new_name)
%             save(fig, fullfile(options.base_path, options.path, append(filename ,  num2str(fignr) , '.m')))
            savefig(fig, fullfile(options.base_path, options.fig_path, append(filename ,  num2str(fignr) , '.fig')))

        elseif options.overwrite
%             save(fullfile(options.base_path, options.path, append(filename , '.m')))
            savefig(fig, fullfile(options.base_path, options.fig_path, append(filename , '.fig')))

        end
    else
        exportgraphics(fig, full_path)
%         save(fullfile(options.base_path, options.path, append(filename, '.m')))
        savefig(fig, fullfile(options.base_path, options.fig_path, append(filename , '.fig')))
    end
end