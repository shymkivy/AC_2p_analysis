function f_save_fig(fig_num, filename, save_path)

% save current fig

if ~numel(fig_num)
    fig = gcf;
else
    fig = figure(fig_num);
end

if ~numel(save_path)
    save_path = 'C:\Users\ys2605\Desktop\stuff\papers\AC_paper';
end

if ~numel(filename)
    num_ch = numel(fig.Children);
    has_title = false(num_ch,1);
    titles_all = cell(num_ch,1);
    for n_ch = 1:num_ch
        if isprop(fig.Children(n_ch), 'Title')
            if isprop(fig.Children(n_ch).Title, 'String')
                if numel(fig.Children(n_ch).Title.String)
                    has_title(n_ch) = 1;
                    titles_all{n_ch} = fig.Children(n_ch).Title.String;
                end
            end
        end
    end

    filename2 = titles_all{find(has_title,1)};
    if iscell(filename2)
        filename2 = cat(2, filename2{:});
    end
    filename2 = regexprep(filename2, ' +', '_');
    filename2 = regexprep(filename2, '[;=.]+', '');

    filename = sprintf('%s_%s_%sh_%sm', filename2, datetime('now', 'Format','yy_MM_dd'), datetime('now', 'Format','hh'), datetime('now', 'Format','mm'));

end


file_path = sprintf('%s\\%s', save_path, filename);

fprintf('saving: %s\n', file_path);

saveas(fig,file_path,'fig')
saveas(fig,file_path,'svg')
saveas(fig,file_path,'png')

fprintf('Done\n')

end





