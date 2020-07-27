function f_plot_cond_reliability(data, ops)



e_colors = {'b', 'r', 'g', 'k'};

for n_tt = 1:numel(ops.dred_params.trial_types_to_dred)
    figure; hold on;
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        cdata = data.(cond_name);
    
        dr_params.cond_name = cond_name;
        dr_params.colors_clust = cat(2,ops.colors_list,ops.colors_list,ops.colors_list);
        dr_params.tt_to_dred_input = ops.dred_params.trial_types_to_dred{n_tt};
        dr_params.num_clust = ops.dred_params.hclust.num_clust{n_tt};
        
        tn_to_dred_all = cell(cdata.num_dsets,1);
        for n_dset = 1:cdata.num_dsets
            [tn_to_dred, trial_type_tag] = f_select_trial_type(dr_params.tt_to_dred_input, cdata, n_dset, ops);
            tn_to_dred_all{n_dset} = repmat(tn_to_dred, cdata.num_cells(n_dset),1);
        end
        tn_to_dred_all = cat(1,tn_to_dred_all{:});
        
        resp_cells = logical(cat(1,cdata.peak_tuned_trials_full{:}));
        cell_reliab = cat(1,cdata.peak_tuned_trials_full_reliab{:});

        trial_num = 30;

        dd_resp_reliab = cell_reliab(resp_cells(:,trial_num),trial_num);

        %[f, x] = ksdensity(dd_resp_reliab, 'Bandwidth',.07);

        %histogram(dd_resp_reliab, 10, 'BinLimits',[0,1], 'FaceColor','none','EdgeColor',e_colors{n_cond},'Normalization','probability', 'LineWidth', 2,'DisplayStyle','stairs')

        [f, x] = ecdf(dd_resp_reliab);
        plot(x,f, 'color', e_colors{n_cond}, 'LineWidth', 2);
    end
    legend(ops.regions_to_analyze, 'Location', 'southeast');
    xlabel('response reliability');
    ylabel('fraction');
    title('Cell reliability distribution');
end


end