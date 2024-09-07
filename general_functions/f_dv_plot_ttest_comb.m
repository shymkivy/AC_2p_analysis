function f_dv_plot_ttest_comb(data, title_tag, ax_lab)

if ~exist('title_tag', 'var') || isempty(title_tag)
    title_tag = '';
end

[num_pts, num_gr] = size(data);

% and t test between dd and cont, and red and cont

h_all = zeros(num_gr, num_gr);
p_all = ones(num_gr, num_gr);
t_all = zeros(num_gr, num_gr);

for ngn1 = 1:num_gr
    for ngn2 = 1:num_gr
        if ngn1>ngn2
            [h,p,~,stats] = ttest(data(:,ngn1), data(:,ngn2), tail='both');
            h_all(ngn1, ngn2) = h;
            p_all(ngn1, ngn2) = p;
            t_all(ngn1, ngn2) = stats.tstat;
        end
    end
end

df = stats.df;

figure; 
sp_all{1} = subplot(1,3,1);imagesc(abs(t_all)); title('t factor');
sp_all{2} = subplot(1,3,2);imagesc(p_all); title('p value, two-tailed');
sp_all{3} = subplot(1,3,3);imagesc(p_all < 0.05); title('is sig');
sgtitle(sprintf('%s;\n paired ttest df=%d', title_tag, df), 'Interpreter', 'none');

if exist('ax_lab', 'var') && ~isempty(ax_lab)
    for n1 = 1:3
        s1 = sp_all{n1};
        s1.XTick = 1:num_gr;
        s1.YTick = 1:num_gr;
        set(s1, 'XTick',s1.XTick, 'XTickLabel',ax_lab) % , 'XTickLabelRotation',30
        set(s1, 'XTick',s1.YTick, 'YTickLabel',ax_lab) % , 'XTickLabelRotation',30
    end
end


%[p,h] = signrank(resp_temp(:,1), resp_temp(:,3))


end