function f_dv_ensless_sintle_trial_corr(app)

n_pl = app.mplSpinner.Value;
tn_list = 3:7;
dset_list = 4:7;

select_resp_cells = app.selectrespcellsCheckBox.Value;
resort_by_ens = app.resortbyensCheckBox.Value;
sort_trials = app.sorttrialsCheckBox.Value;
sort_with_full_firing_rate = app.sortwithfullfrCheckBox.Value;

corr_vals = zeros(numel(dset_list), numel(tn_list));
for n_dset = 1:numel(dset_list)
    for n_tn = 1:numel(tn_list)
        dset_idx = dset_list(n_dset);
        
        ddata = app.data(dset_idx,:);
        tn_all = tn_list(n_tn);
        
        if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
            cdata = ddata.cdata{n_pl};
        else
            cdata = cat(1,ddata.cdata{:});
        end
        
        num_cells = sum([cdata.num_cells]);
        firing_rate = cat(1,cdata.S_sm);
        
        stats1 = ddata.stats;
        if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
            resp_cell = logical(sum(stats1{n_pl}.resp_cells_peak(:,tn_all),2));
        else
            resp_cells_all = cell(numel(stats1),1);
            for n_pl2 = 1:numel(stats1)
                resp_cells_all{n_pl2} = logical(sum(stats1{n_pl2}.resp_cells_peak(:,tn_all),2));
            end
            resp_cell = cat(1,resp_cells_all{:});
        end
        
        %%
        trial_types = ddata.trial_types{1};
        stim_frame_index = ddata.stim_frame_index{1};
        
        trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
        [~, trial_frames] = f_dv_compute_window_t(trial_window, cdata.volume_period);

        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_frames);
        if ~isempty(ddata.MMN_freq{1})
            [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, ddata.MMN_freq{1}, app.ops);
        else
            trial_data_sort_wctx = trial_data_sort;
            trial_types_wctx = trial_types;
        end

        tr_idx = logical(sum(trial_types_wctx == app.ops.context_types_all(tn_all)',2));
        tr_loc_full = find(tr_idx);
        tr_data = trial_data_sort_wctx(:,:,tr_idx);

        if select_resp_cells
            sel_resp_cells = logical(sum(resp_cell,2));
            tr_data2 = tr_data(sel_resp_cells,:,:);
            resp_all2 = resp_cell(sel_resp_cells,:);
            firing_rate2 = firing_rate(sel_resp_cells,:);
        else
            tr_data2 = tr_data;
            resp_all2 = resp_cell;
            firing_rate2 = firing_rate;
        end

        [num_cells2, ~, num_tr] = size(tr_data2);

        hc_params.plot_dist_mat = 0;
        hc_params.plot_clusters = 0;
        hc_params.num_clust = 1;


        tr_data_2d_tr = reshape(tr_data2, [], num_tr);
        hclust_out_trial = f_hcluster_wrap(tr_data_2d_tr', hc_params);
        tr_data3 = tr_data2(:,:,hclust_out_trial.dend_order);

        SI = 1-hclust_out_trial.dist;
        SI_vals = tril(SI,-1);
        SI_vals(SI_vals==0) = [];
        
        corr_vals(n_dset, n_tn) = mean(SI_vals);
   
    end
end

color1 = jet(10);
figure; hold on
for n_freq = 1:numel(tn_list)
    plot([0.5 1 2 4], corr_vals(:,n_freq), 'o-', 'color', color1(tn_list(n_freq),:), 'linewidth', 2)
end
xlabel('ISI duration'); ylabel('Pairwise correlation')
title('Mean pairwise correlations');

figure; imagesc(reshape(color1, [10 1 3]))

end