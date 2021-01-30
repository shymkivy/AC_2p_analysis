function f_mpl_trial_trial_dist(data, ops)

plot_stuff = 0;
trial_mean_vec = cell(numel(ops.regions_to_analyze),1);
trial_raster = cell(numel(ops.regions_to_analyze),1);
tt_tag = cell(numel(ops.dred_params.trial_types_for_dist),1);
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    num_dsets = numel(cdata.area);
    
    dr_params.cond_name = cond_name;
    trial_mean_vec_tt = cell(numel(ops.dred_params.trial_types_for_dist),num_dsets);
    trial_raster_tt = cell(numel(ops.dred_params.trial_types_for_dist),num_dsets);
    for n_tt = 1:numel(ops.dred_params.trial_types_for_dist)
        
        for n_dset = 1:num_dsets
            %%
            dr_params.tt_to_dred_input = ops.dred_params.trial_types_for_dist{n_tt};
            disp([cond_name, ' dset ' num2str(n_dset)]);
            trial_types = cdata.trial_types_wctx{n_dset};
            trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
            trial_data_sort_sm_wctx = cdata.trial_data_sort_sm_wctx{n_dset};

            %% select trials
            [tn_to_dred, trial_type_tag] = f_select_trial_type(dr_params.tt_to_dred_input, cdata, n_dset, ops);
            tt_to_dred = ops.context_types_all(tn_to_dred);
            tt_tag{n_tt} = trial_type_tag;
            %%
            trials_idx_dred = logical(sum(trial_types == tt_to_dred(:)' ,2));
            trial_peaks_dred = trial_peaks(:,trials_idx_dred);
            trial_types_dred = trial_types(trials_idx_dred);
            trial_data_sort_sm_wctx = trial_data_sort_sm_wctx(:,:,trials_idx_dred);
            
            extra_red = logical((trial_types==160)+(trial_types==260));
            trial_peaks_cut = trial_peaks(:, ~extra_red);
            trial_types_cut = trial_types(~extra_red);
            
            if ops.dred_params.use_responsive_cells
                resp_cells = cdata.peak_tuned_trials_combined{n_dset};
                %resp_cells = logical(sum(resp_cells(:,tn_to_dred),2));
                resp_cells = logical(sum(resp_cells,2));
                trial_peaks_dred = trial_peaks_dred(resp_cells,:);
                trial_peaks_cut = trial_peaks_cut(resp_cells,:);
                trial_data_sort_sm_wctx = trial_data_sort_sm_wctx(resp_cells,:,:);
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
            [~, trial_ind_cut] = sort(trial_types_cut);
            
            %%
            trial_peaks_dred_sort_cell = trial_peaks_dred(hclust_out_cell.dend_order,:);
            trial_peaks_dred_sort = trial_peaks_dred_sort_cell(:,trial_ind);
            trial_types_dred_sort = trial_types_dred(trial_ind);
            
            trial_peaks_cut_sort_cell = trial_peaks_cut(hclust_out_cell.dend_order,:);
            trial_peaks_cut_sort = trial_peaks_cut_sort_cell(:,trial_ind_cut);
            trial_types_cut_sort = trial_types_cut(trial_ind_cut);
            
            if plot_stuff
                
                rand_seq = randsample(num_cells, num_cells);
                trial_num_cut = f_tt_to_tn(trial_types_cut, ops,1);
                f_plot_raster_mean(trial_peaks_cut(rand_seq,:), 1,trial_num_cut, ops)
                title('raster');

                f_plot_raster_mean(trial_peaks_cut_sort_cell, 1,trial_num_cut, ops)
                title('raster cell sorted');
                
                trial_num_cut_sort = f_tt_to_tn(trial_types_cut_sort, ops,1);
                f_plot_raster_mean(trial_peaks_cut_sort,trial_num_cut_sort, ops, 1)
                title('raster cell trial sorted');
                
                trial_num_dred_sort = f_tt_to_tn(trial_types_dred_sort, ops,1);
                f_plot_raster_mean(trial_peaks_dred_sort, 1,trial_num_dred_sort, ops)
                title('raster cell freq trial sorted');
                
            end
            
            
            trial_mean_vec_tt{n_tt,n_dset} = zeros(num_cells, numel(tt_to_dred));
            trial_raster_tt{n_tt,n_dset} = cell(numel(tt_to_dred),1);
            
            for n_tr = 1:numel(tt_to_dred)
                temp_ras = trial_peaks_dred_sort(:,trial_types_dred_sort == tt_to_dred(n_tr));
                trial_raster_tt{n_tt,n_dset}{n_tr} = temp_ras;
                temp_mean_vec = mean(temp_ras,2);
                trial_mean_vec_tt{n_tt,n_dset}(:,n_tr) = temp_mean_vec;
                
            end
            
            if plot_stuff
                scale_fac = .9/max(trial_mean_vec_tt{n_tt,n_dset}(:));
                f1 = figure; hold on; axis tight;
                for n_tr = 1:numel(tt_to_dred)
                    plot(trial_mean_vec_tt{n_tt,n_dset}(:,n_tr)*scale_fac+n_tr, 1:num_cells, 'LineWidth', 2, 'color', ops.context_types_all_colors(n_tr,:));
                    line([n_tr n_tr], [1 num_cells], 'color', 'k');
                    
                end
                title(sprintf('%s, %s, dset%d population vec',trial_type_tag,  cond_name, n_dset));
                f1.Children.YAxis.Direction = 'reverse';
            end
            
%             figure; imagesc(trial_ave_mat)
%             title(sprintf('%s, %s, dset%d',trial_type_tag,  cond_name, n_dset));
            
            if plot_stuff
                corr_sim = trial_mean_vec_tt{n_tt,n_dset}' * trial_mean_vec_tt{n_tt,n_dset};
                figure; imagesc(corr_sim);
                title(sprintf('%s, %s, dset%d, prod',trial_type_tag,  cond_name, n_dset));
            end
            if plot_stuff
                % same as cosine with normalized traces
                corr_sim = corr(trial_mean_vec_tt{n_tt,n_dset},trial_mean_vec_tt{n_tt,n_dset});
                figure; imagesc(corr_sim);
                title(sprintf('%s, %s, dset%d, corr',trial_type_tag,  cond_name, n_dset));
            end

        end
    end
    trial_mean_vec{n_cond} = trial_mean_vec_tt;
    trial_raster{n_cond} =  trial_raster_tt;
end

%% plot all
f_plot_stim_vec_dist_mat(trial_mean_vec, tt_tag, ops)

%% mmn plots
if sum(strcmpi(ops.dred_params.trial_types_for_dist, 'mmn12'))
    tt_ind = find(strcmpi(ops.dred_params.trial_types_for_dist, 'mmn12'));
    sim_ind = [1,3; 4,6];
    f_plot_stim_vec_dist(trial_mean_vec, tt_ind, sim_ind, ops, 'cosine'); %cosineSI, cosine, euclidean
    title('cont-dd')
    %ylim([0 .6]);
    
    f_plot_stim_vec_trial_dist(trial_mean_vec, trial_raster, ops, 'cosine'); %cosineSI, cosine, euclidean
    

    sim_ind = [2,3; 5,6];
    f_plot_stim_vec_dist(trial_mean_vec, tt_ind, sim_ind, ops, 'cosine');
    title('red-dd')

    sim_ind = [1,3; 4,6];
    f_plot_stim_vec_dist_v_ncell(trial_mean_vec, tt_ind, sim_ind, ops, 99999);
    title('cont-dd')
end

%% freq plots
tt_ind = 1;
sim_ind = [3,4; 4,5; 5,6; 6,7; 7,8];
%sim_ind = [1,3; 2,4; 3,5; 4,6; 5,7; 6,8; 7,9; 8,10];

f_plot_stim_vec_dist(trial_mean_vec, tt_ind, sim_ind, ops);
title('adj freqs')

f_plot_stim_vec_dist(trial_mean_vec, tt_ind, sim_ind, ops);
title('adj freqs 2 apart')




end