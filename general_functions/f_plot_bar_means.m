function f_plot_bar_means(means, sems, labels, colors_big, colors_small)

[num_sbar, num_bbar] = size(means);

figure; hold on;
bar(categorical(labels(:,1), labels(:,1)), [0 0 0]);
for n_bb = 1:num_bbar
    b1 = bar(n_bb, means(:,n_bb));
    num_reg2 = numel(means(:,n_bb));
    if num_reg2 > 1
        errorbar(n_bb + ((1:4)-2.5)/5.5, means(:,n_bb), sems(:,n_bb), '.k');
        for n_sb = 1:num_sbar
            b1(n_sb).FaceColor = colors_small{n_sb};
        end
    else
        errorbar(n_bb, means(:,n_bb), sems(:,n_bb), '.k');
        b1.FaceColor = colors_big{n_bb};
    end
end


end