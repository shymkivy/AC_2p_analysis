function p_vals_adj = f_FDR_correction(p_vals)

[~, s_idx] = sort(p_vals);
 [~, s_idx_rev] = sort(s_idx);

num_vals = numel(p_vals);
rank1 = (1:num_vals)';
rank2 = rank1(s_idx_rev);

p_vals4 = p_vals*num_vals./rank2;

p_vals_adj = p_vals4;
last_val = 1;
for n_val = 1:num_vals
    idx1 = rank2 == (num_vals-n_val+1);
    p_vals_adj(idx1) = min([p_vals_adj(idx1), last_val]);
    last_val = p_vals_adj(idx1);
end

end