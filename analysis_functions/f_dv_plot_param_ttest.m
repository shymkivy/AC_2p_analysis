function f_dv_plot_param_ttest(means, stds, dfs, title_tag, ax_lab)
% fisher lsd t = (mean1 - mean2) / sqrt(MSE(1/n1 - 1/n2))

% t = (tau1 - tau2)/sqrt(var1 + var2)
% df = df1 + df2

diff2 = abs(means - means');

stds2 = sqrt(stds.^2 + stds'.^2);

t_vals1 = tril(diff2./stds2);

p_vals1 = (1 - tcdf(t_vals1, dfs+dfs'))*2;

[d1, d2] = size(t_vals1);

sp_all = cell(3,1);

figure; 
sp_all{1} = subplot(1,3,1);imagesc(abs(t_vals1)); title('t factor');
sp_all{2} = subplot(1,3,2);imagesc(p_vals1); title('p value, two-tailed');
sp_all{3} = subplot(1,3,3);imagesc(p_vals1 < 0.05); title('is sig');
sgtitle(sprintf('%s; ttest', title_tag), 'Interpreter', 'none');

if exist('ax_lab', 'var') && ~isempty(ax_lab)
    for n1 = 1:3
        s1 = sp_all{n1};
        s1.XTick = 1:d2;
        s1.YTick = 1:d1;
        set(s1, 'XTick',s1.XTick, 'XTickLabel',ax_lab) % , 'XTickLabelRotation',30
        set(s1, 'XTick',s1.YTick, 'YTickLabel',ax_lab) % , 'XTickLabelRotation',30
    end
end


end