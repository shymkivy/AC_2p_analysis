function f_dv_plot_anova1(p_val, tbl, stats, title_tag, ax_lab, paired, do_fdr)
% fisher lsd t = (mean1 - mean2) / sqrt(MSE(1/n1 - 1/n2))

if ~exist('do_fdr', 'var')
    do_fdr = 1;
end

if ~exist('paired', 'var')
    paired = 0;
end

if numel(unique(stats.n)) > 1
    paired = 0;
end

if paired
    t_tag = 'paired';
else
    t_tag = ''; 
end

if (p_val < 0.05) && ~isnan(stats.s)

    idx1 = strcmpi(tbl(1,:), 'MS');
    idx2 = strcmpi(tbl(1,:), 'df');
    idx3 = strcmpi(tbl(1,:), 'F');
    idx4 = or(strcmpi(tbl(:,1), 'Groups'), strcmpi(tbl(:,1), 'Columns'));
    idx5 = strcmpi(tbl(:,1), 'Error');
    
    MSE = tbl{idx5,idx1};
    dfc = tbl{idx4,idx2};
    dfe = tbl{idx5,idx2};
    Fc = tbl{idx4,idx3};
    
    mean_diff = stats.means - stats.means';
    [d1, d2] = size(mean_diff);

    if paired
        df_cat = (stats.n + stats.n')/2-1;
        den_t = MSE./df_cat;
    else
        df_cat = stats.n + stats.n' - 1;
        den_t = MSE./stats.n + MSE./stats.n';
    end

    %t_vals1 = tril(mean_diff)./sqrt(MSE.*(1./(stats.n-1) + 1./(stats.n-1)'));
    t_vals1 = tril(mean_diff)./sqrt(den_t);

    p_vals1 = (1 - tcdf(abs(t_vals1), df_cat))*2;
    
    p_vals2 = tril(p_vals1, -1);

    p_vals3 = p_vals2(p_vals2>0);
    
    p_vals_adj2 = f_FDR_correction(p_vals3);

    p_vals_adj = p_vals1;
    p_vals_adj(p_vals2>0) = p_vals_adj2;

    sp_all = cell(5,1);

    figure; 
    sp_all{1} = subplot(2,3,1);imagesc(abs(t_vals1)); title('t factor');
    sp_all{2} = subplot(2,3,2);imagesc(p_vals1); title('p value, two-tailed');
    sp_all{3} = subplot(2,3,3);imagesc(p_vals1 < 0.05); title('is sig');
    sp_all{4} = subplot(2,3,4);imagesc(df_cat); title('df');
    sp_all{5} = subplot(2,3,5);imagesc(p_vals_adj); title('p value FDR adj, two-tailed');
    sp_all{6} = subplot(2,3,6);imagesc(p_vals_adj < 0.05); title('is sig FDR adj');
    sgtitle(sprintf('%s;\n anova p=%.4g; Fc=%.4g; dfe=%d; dfc=%d; t%s', title_tag, p_val, Fc, dfe, dfc, t_tag), 'Interpreter', 'none');
    
    if ~exist('ax_lab', 'var') || isempty(ax_lab)
        ax_lab = cellstr(num2str((1:numel(stats.means))'));
    end
    for n1 = 1:6
        s1 = sp_all{n1};
        s1.XTick = 1:d2;
        s1.YTick = 1:d1;
        set(s1, 'XTick',s1.XTick, 'XTickLabel',ax_lab) % , 'XTickLabelRotation',30
        set(s1, 'XTick',s1.YTick, 'YTickLabel',ax_lab) % , 'XTickLabelRotation',30
    end
    
    fprintf('%s\n', title_tag);
    if p_val>0.01
        sig_tag = '*';
    elseif p_val>0.001
        sig_tag = '**';
    elseif p_val<=0.001
        sig_tag = '***';
    end
    fprintf('Fanova(%d)=%.2f, p=%.2e%s\n', dfe, Fc, p_val, sig_tag)
    if do_fdr
        if_print_p(t_vals1, p_vals_adj, stats, df_cat, ax_lab, 1, t_tag)
    else
        if_print_p(t_vals1, p_vals1, stats, df_cat, ax_lab, 0, t_tag)
    end

end

end

function if_print_p(t_vals1, p_vals1, stats, df_cat, ax_lab, is_adj, t_tag)

t_vals_flat = abs(t_vals1(:));

[d1, d2] = size(t_vals1);

[t_sort, idx1] = sort(t_vals_flat, 'descend');
vals_flat = 1:(d1*d2);
vals_sort = vals_flat(idx1);
p_flat = p_vals1(:);
p_sort = p_flat(idx1);
df_cat_sort = df_cat(idx1);

for n_st = 1:numel(t_sort)
    if p_sort(n_st) < 0.05
        [x1, y1] = ind2sub([d1, d2], vals_sort(n_st));
        if mean(stats.means) < 0.1
            if x1 < y1
                if stats.means(x1) > stats.means(y1)
                    tag2 = '>';
                else
                    tag2 = '<';
                end
                tag1 = sprintf('%s(%.2e) %s %s(%.2e)', ax_lab{x1}, stats.means(x1), tag2, ax_lab{y1}, stats.means(y1));
            else
                if stats.means(y1) > stats.means(x1)
                    tag2 = '>';
                else
                    tag2 = '<';
                end
                tag1 = sprintf('%s(%.2e) %s %s(%.2e)', ax_lab{y1}, stats.means(y1), tag2, ax_lab{x1}, stats.means(x1));
            end
        else
            if x1 < y1
                if stats.means(x1) > stats.means(y1)
                    tag2 = '>';
                else
                    tag2 = '<';
                end
                tag1 = sprintf('%s(%.2f) %s %s(%.2f)', ax_lab{x1}, stats.means(x1), tag2, ax_lab{y1}, stats.means(y1));
            else
                if stats.means(y1) > stats.means(x1)
                    tag2 = '>';
                else
                    tag2 = '<';
                end
                tag1 = sprintf('%s(%.2f) %s %s(%.2f)', ax_lab{y1}, stats.means(y1), tag2, ax_lab{x1}, stats.means(x1));
            end
        end
        if p_sort(n_st)>0.01
            sig_tag = '*';
        elseif p_sort(n_st)>0.001
            sig_tag = '**';
        elseif p_sort(n_st)<=0.001
            sig_tag = '***';
        end

        if is_adj
            fprintf('%s, t%s(%d)=%.2f, p_adj=%.2e%s\n', tag1, t_tag, df_cat_sort(n_st), t_sort(n_st), p_sort(n_st), sig_tag)
        else
            fprintf('%s, t%s(%d)=%.2f, p=%.2e%s\n', tag1, t_tag, df_cat_sort(n_st), t_sort(n_st), p_sort(n_st), sig_tag)
        end
    end
end
fprintf('\n');

end
