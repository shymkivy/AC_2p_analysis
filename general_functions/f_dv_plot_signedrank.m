function f_dv_plot_signedrank(data, title_tag, ax_lab)

% nonparametric equivalent of paired ttest

if ~exist('title_tag', 'var') || isempty(title_tag)
    title_tag = '';
end

[num_pts, num_gr] = size(data);

% and t test between dd and cont, and red and cont

h_all = zeros(num_gr, num_gr);
p_all = ones(num_gr, num_gr);
z_all = zeros(num_gr, num_gr);

for ngn1 = 1:num_gr
    for ngn2 = 1:num_gr
        if ngn1>ngn2
            [p,h,stats] = signrank(data(:,ngn1), data(:,ngn2), tail='both');
            h_all(ngn1, ngn2) = h;
            p_all(ngn1, ngn2) = p;
            z_all(ngn1, ngn2) = stats.zval;
        end
    end
end

df = num_pts-1;

figure; 
sp_all{1} = subplot(1,3,1);imagesc(abs(z_all)); title('z value');
sp_all{2} = subplot(1,3,2);imagesc(p_all); title('p value, two-tailed');
sp_all{3} = subplot(1,3,3);imagesc(p_all < 0.05); title('is sig');
sgtitle(sprintf('%s;\n signed-rank test df=%d', title_tag, df), 'Interpreter', 'none');

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