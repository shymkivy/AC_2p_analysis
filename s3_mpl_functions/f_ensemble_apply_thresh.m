function ens_out = f_ensemble_apply_thresh(coeffs, scores, thresh_coeffs, thresh_scores, num_ens, params)
plot_stuff = f_get_param(params, 'plot_stuff', 0);

%% first detect cells

two_sided = logical(size(thresh_coeffs,2)-1)

for n_ens = 1:num_ens
    factors1 = coeffs(:,n_ens);
    [fac_sort, indx1] = sort(factors1, 'ascend');
    
    ens_cells1 = indx1(fac_sort>max(thresh_coeffs(n_ens,:)))
    if two_sided
        ens_cells2 = indx1(fac_sort<min(thresh_coeffs(n_ens,:)))
    end
end

X = coeffs';
comps_cell = if_get_comp_thresh(X, num_ens, 2, plot_stuff);
if plot_stuff
    title('cells');
end

X = scores;
comps_tr = if_get_comp_thresh(X, num_ens, 2, plot_stuff);
if plot_stuff
    title('trials');
end

ens_out.trials = comps_tr;
ens_out.cells = comps_cell;

end
