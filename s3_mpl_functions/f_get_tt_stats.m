function p_val_mat = f_get_tt_stats(corr_list)

num_groups = numel(corr_list);
p_val_mat = zeros(num_groups,num_groups);

for n_gr1 = 1:num_groups
    for n_gr2 = 1:num_groups
        [~, p] = ttest2(corr_list{n_gr1}, corr_list{n_gr2});
        p_val_mat(n_gr1, n_gr2) = p;
    end
end

end