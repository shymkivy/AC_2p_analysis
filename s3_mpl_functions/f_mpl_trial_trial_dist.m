function f_mpl_trial_trial_dist(data, ops)

plot_stuff = 0;
trial_ave_vec = cell(numel(ops.regions_to_analyze),1);
tt_tag = cell(numel(ops.dred_params.trial_types_for_dist),1);
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    dr_params.cond_name = cond_name;
    trial_ave_vec_tt = cell(numel(ops.dred_params.trial_types_for_dist),numel(cdata.num_dsets,1));
    for n_tt = 1:numel(ops.dred_params.trial_types_for_dist)
        
        for n_dset = 1:cdata.num_dsets
            %%
            dr_params.tt_to_dred_input = ops.dred_params.trial_types_for_dist{n_tt};
            disp([cond_name, ' dset ' num2str(n_dset)]);
            trial_types = cdata.trial_types_pr{n_dset};
            trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
            trial_data_sort_sm_pr = cdata.trial_data_sort_sm_pr{n_dset};

            %% select trials
            [tn_to_dred, trial_type_tag] = f_select_trial_type(dr_params.tt_to_dred_input, cdata, n_dset, ops);
            tt_to_dred = ops.context_types_all(tn_to_dred);
            tt_tag{n_tt} = trial_type_tag;
            %%
            trials_idx_dred = logical(sum(trial_types == tt_to_dred(:)' ,2));
            trial_peaks_dred = trial_peaks(:,trials_idx_dred);
            trial_types_dred = trial_types(trials_idx_dred);
            trial_data_sort_sm_pr = trial_data_sort_sm_pr(:,:,trials_idx_dred);

            if ops.dred_params.use_responsive_cells
                resp_cells = cdata.peak_tuned_trials_combined{n_dset};
                resp_cells = logical(sum(resp_cells(:,tn_to_dred),2));
                trial_peaks_dred = trial_peaks_dred(resp_cells,:);
                trial_data_sort_sm_pr = trial_data_sort_sm_pr(resp_cells,:,:);
            end
            
            num_cells = size(trial_peaks_dred,1);

            %% cluster cell types for aesthetics
            hc_params.method = ops.dred_params.hclust.method;
            hc_params.metric = ops.dred_params.hclust.plot_metric;
            hc_params.plot_dist_mat = 0;
            hc_params.plot_clusters = 0;
            hclust_out_cell = f_hcluster_cell(trial_peaks_dred, [], hc_params, ops);
            
            %% sort trial types
            [~, trial_ind] = sort(trial_types_dred);
            
            %%
            trial_peaks_dred_sort = trial_peaks_dred(hclust_out_cell.dend_order,trial_ind);
            trial_types_dred_sort = trial_types_dred(trial_ind);
            if plot_stuff
                figure; imagesc(trial_types_dred_sort, 1:num_cells, trial_peaks_dred_sort)
            end
            
            trial_ave_vec_tt{n_tt,n_dset} = zeros(num_cells, numel(tt_to_dred));
            
            for n_tr = 1:numel(tt_to_dred)
                temp_ras = trial_peaks_dred_sort(:,trial_types_dred_sort == tt_to_dred(n_tr));
                trial_ave_vec_tt{n_tt,n_dset}(:,n_tr) = mean(temp_ras,2);
            end
            scale_fac = .9/max(trial_ave_vec_tt{n_tt,n_dset}(:));
            if plot_stuff
                figure; hold on; axis tight;
                for n_tr = 1:numel(tt_to_dred)
                    plot(trial_ave_vec_tt{n_tt,n_dset}(:,n_tr)*scale_fac+n_tr, 1:num_cells, 'LineWidth', 2);
                    line([n_tr n_tr], [1 num_cells], 'color', 'k')
                end
                title(sprintf('%s, %s, dset%d population vec',trial_type_tag,  cond_name, n_dset));
            end
%             figure; imagesc(trial_ave_mat)
%             title(sprintf('%s, %s, dset%d',trial_type_tag,  cond_name, n_dset));
            
            if plot_stuff
                corr_sim = trial_ave_vec_tt{n_tt,n_dset}' * trial_ave_vec_tt{n_tt,n_dset};
                figure; imagesc(corr_sim);
                title(sprintf('%s, %s, dset%d, prod',trial_type_tag,  cond_name, n_dset));
            end
            if plot_stuff
                % same as cosine with normalized traces
                corr_sim = corr(trial_ave_vec_tt{n_tt,n_dset},trial_ave_vec_tt{n_tt,n_dset});
                figure; imagesc(corr_sim);
                title(sprintf('%s, %s, dset%d, corr',trial_type_tag,  cond_name, n_dset));
            end

        end
    end
    trial_ave_vec{n_cond} = trial_ave_vec_tt;
end

%% plot all
f_plot_corr_mat(trial_ave_vec, tt_tag, ops)

tt_ind = find(strcmpi(ops.dred_params.trial_types_for_dist, 'mmn12'));
sim_ind = [1,3; 4,6];
f_plot_corr_stim(trial_ave_vec, tt_ind, sim_ind, ops);
title('cont-dd')

tt_ind = 1;
sim_ind = [1,2; 2,3; 3,4; 4,5; 5,6; 6,7; 7,8; 8,9; 9,10];
f_plot_corr_stim(trial_ave_vec, tt_ind, sim_ind, ops);
title('adj freqs')

tt_ind = 1;
sim_ind = [1,3; 2,4; 3,5; 4,6; 5,7; 6,8; 7,9; 8,10;];
f_plot_corr_stim(trial_ave_vec, tt_ind, sim_ind, ops);
title('adj freqs')


end