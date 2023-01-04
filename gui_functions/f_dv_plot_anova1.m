function f_dv_plot_anova1(p_val, tbl, stats, title_tag)
% fisher lsd t = (mean1 - mean2) / sqrt(MSE(1/n1 - 1/n2))

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

    t_vals1 = tril(stats.means - stats.means')./sqrt(MSE.*(1./stats.n + 1./stats.n'));

    p_vals1 = (1 - tcdf(abs(t_vals1), dfe))*2;

    figure; 
    subplot(1,3,1);imagesc(abs(t_vals1)); title('t factor');
    subplot(1,3,2);imagesc(p_vals1); title('p value, two-tailed');
    subplot(1,3,3);imagesc(p_vals1 < 0.05); title('is sig');
    sgtitle(sprintf('%s;\n anova p=%.4g; Fc=%.4g; dfe=%d; dfc=%d', title_tag, p_val, Fc, dfe, dfc), 'Interpreter', 'none')
end

end