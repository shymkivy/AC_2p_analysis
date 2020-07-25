function f_mpl_population_analysis_trials3(data, ops)

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    %% cycle through analysis groups 
    for n_tt = 1:numel(ops.dred_params.trial_types_to_dred)

        dr_params.cond_name = cond_name;
        dr_params.colors_clust = cat(2,ops.colors_list,ops.colors_list,ops.colors_list);
        dr_params.tt_to_dred_input = ops.dred_params.trial_types_to_dred{n_tt};
        dr_params.num_clust = ops.dred_params.hclust.num_clust{n_tt};
        
        if ops.dred_params.do_hclust
            f_hcluster_cond(cdata, dr_params, ops);      
        end
    end
end

end

