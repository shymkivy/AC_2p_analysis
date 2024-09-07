function f_plot_bar(means, labels, colors_big, colors_small)

[num_sbar, num_bbar] = size(means);

figure; hold on;
bar(categorical(labels(:,1), labels(:,1)), zeros(num_bbar,1));

means1 = zeros(num_sbar, num_bbar);
sems1 = zeros(num_sbar, num_bbar);
for n_bb = 1:num_bbar
    num_pts = size(means{1,n_bb},1);
    for n_sb = 1:num_sbar
        temp_data = means{n_sb,n_bb};
        means1(n_sb,n_bb) = mean(temp_data,1);
        sems1(n_sb,n_bb) = std(temp_data,[],1)/sqrt(num_pts-1);
    end
end
    
for n_bb = 1:num_bbar
    b1 = bar(n_bb, means1(:,n_bb));
    if num_sbar > 1
        errorbar(n_bb + ((1:4)-2.5)/5.5, means1(:,n_bb), sems1(:,n_bb), '.k');
        for n_sb = 1:num_sbar
            b1(n_sb).FaceColor = colors_small{n_sb};
        end
    else
        errorbar(n_bb, means1(:,n_bb), sems1(:,n_bb), '.k');
        b1.FaceColor = colors_big{n_bb};
    end
end


end