function f_dv_plot_anovan(p_val, tbl, stats, data, lab_gr, title_tag)
% fisher lsd t = (mean1 - mean2) / sqrt(MSE(1/n1 - 1/n2))




if sum(p_val < 0.05)
    
    num_gr = numel(lab_gr);
    means_gr = cell(num_gr,1);

    for n_gr = 1:num_gr
        gr1 = lab_gr{n_gr};
        gr_lab = unique(gr1);
        num_lab = numel(gr_lab);
        means2 = zeros(num_lab, 1);
        for n_lab = 1:num_lab
            idx1 = strcmpi(gr1, gr_lab{n_lab});
            means2(n_lab) = mean(data(idx1));
        end
        means_gr{n_gr} = means2;
    end

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

    t_vals1 = tril(mean_diff)./sqrt(MSE.*(1./(stats.n-1) + 1./(stats.n-1)'));

    p_vals1 = (1 - tcdf(abs(t_vals1), dfe))*2;
    
    [d1, d2] = size(t_vals1);
    
    sp_all = cell(3,1);

    figure; 
    sp_all{1} = subplot(1,3,1);imagesc(abs(t_vals1)); title('t factor');
    sp_all{2} = subplot(1,3,2);imagesc(p_vals1); title('p value, two-tailed');
    sp_all{3} = subplot(1,3,3);imagesc(p_vals1 < 0.05); title('is sig');
    sgtitle(sprintf('%s;\n anova p=%.4g; Fc=%.4g; dfe=%d; dfc=%d', title_tag, p_val, Fc, dfe, dfc), 'Interpreter', 'none');
    
    if ~exist('ax_lab', 'var') || isempty(ax_lab)
        ax_lab = cellstr(num2str((1:numel(stats.means))'));
    end
    for n1 = 1:3
        s1 = sp_all{n1};
        s1.XTick = 1:d2;
        s1.YTick = 1:d1;
        set(s1, 'XTick',s1.XTick, 'XTickLabel',ax_lab) % , 'XTickLabelRotation',30
        set(s1, 'XTick',s1.YTick, 'YTickLabel',ax_lab) % , 'XTickLabelRotation',30
    end
    
    t_vals_flat = abs(t_vals1(:));
    
    [t_sort, idx1] = sort(t_vals_flat, 'descend');
    vals_flat = 1:(d1*d2);
    vals_sort = vals_flat(idx1);
    p_flat = p_vals1(:);
    p_sort = p_flat(idx1);
    
    fprintf('%s\n', title_tag);
    if p_val>0.01
        sig_tag = '*';
    elseif p_val>0.001
        sig_tag = '**';
    elseif p_val<=0.001
        sig_tag = '***';
    end
    fprintf('Fanova(%d)=%.2f, p=%.2e%s\n', dfe, Fc, p_val, sig_tag)
 
    for n_st = 1:numel(t_sort)
        if p_sort(n_st) < 0.05
            [x1, y1] = ind2sub([d1, d2], vals_sort(n_st));
            if mean(stats.means) < 0.1
                if x1 < y1
                    if stats.means(x1) > stats.means(y1)
                        tag1 = sprintf('%s(%.2e)>%s(%.2e)', ax_lab{x1}, stats.means(x1), ax_lab{y1}, stats.means(y1));
                    else
                        tag1 = sprintf('%s(%.2e)<%s(%.2e)', ax_lab{x1}, stats.means(x1), ax_lab{y1}, stats.means(y1));
                    end
                else
                    if stats.means(y1) > stats.means(x1)
                        tag1 = sprintf('%s(%.2e)>%s(%.2e)', ax_lab{y1}, stats.means(y1), ax_lab{x1}, stats.means(x1));
                    else
                        tag1 = sprintf('%s(%.2e)<%s(%.2e)', ax_lab{y1}, stats.means(y1), ax_lab{x1}, stats.means(x1));
                    end
                end
            else
                if x1 < y1
                    if stats.means(x1) > stats.means(y1)
                        tag1 = sprintf('%s(%.2f)>%s(%.2f)', ax_lab{x1}, stats.means(x1), ax_lab{y1}, stats.means(y1));
                    else
                        tag1 = sprintf('%s(%.2f)<%s(%.2f)', ax_lab{x1}, stats.means(x1), ax_lab{y1}, stats.means(y1));
                    end
                else
                    if stats.means(y1) > stats.means(x1)
                        tag1 = sprintf('%s(%.2f)>%s(%.2f)', ax_lab{y1}, stats.means(y1), ax_lab{x1}, stats.means(x1));
                    else
                        tag1 = sprintf('%s(%.2f)<%s(%.2f)', ax_lab{y1}, stats.means(y1), ax_lab{x1}, stats.means(x1));
                    end
                end
            end
            if p_sort(n_st)>0.01
                sig_tag = '*';
            elseif p_sort(n_st)>0.001
                sig_tag = '**';
            elseif p_sort(n_st)<=0.001
                sig_tag = '***';
            end


            fprintf('%s, tpaired(%d)=%.2f, p=%.2e%s\n', tag1, dfe, t_sort(n_st), p_sort(n_st), sig_tag)
        end
    end
    fprintf('\n');

end

end

% 
% df = 2;
% rest_tempn = resp_temp-mean(resp_temp,1);
% 
% SSR = sum((mean(resp_temp,1) - mean(mean(resp_temp,1))).^2)*2837;
% 
% SSE = sum(rest_tempn(:).^2);
% 
% SST = var(resp_temp(:))*2837*3;
% 
% MSE = SSE/(2837*3-df-1);
% 
% MSR = SSR/df;
% 
% MSR/MSE
