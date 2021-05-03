function f_mpl_population_analysis_trials3(data, ops)

data_out = struct();
dr_params.dim_est_st = struct('cond_name', [], 'n_dset', [], 'num_cells', [],...
    'num_cells_samp', [], 'num_comp_est', [], 'n_rep', []);
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    
    %% cycle through analysis groups 
    for n_tt = 1:numel(ops.dred_params.trial_types_to_dred)
        
        dr_params.cond_name = cond_name;
        dr_params.colors_clust = cat(2,ops.colors_list,ops.colors_list,ops.colors_list,ops.colors_list,ops.colors_list);
        dr_params.tt_to_dred_input = ops.dred_params.trial_types_to_dred{n_tt};
        dr_params.num_clust = ops.dred_params.hclust.num_clust{n_tt};
        
        if ops.dred_params.do_hclust
            [dr_params, hdata_out] = f_hcluster_cond(cdata, dr_params, ops);      
        end
        data_out.(dr_params.trial_type_tag).(cond_name) = hdata_out;
    end
end

f_plot_mean_similarity(data_out, ops);

f_plot_dset_dimensionality(data_out, ops);

if ops.dred_params.do_dim_estimate
     f_plot_dim_est_data(dr_params.dim_est_st, ops)
end

end

