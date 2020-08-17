function [dr_params, hdata_out] = f_hcluster_cond(cdata, dr_params, ops)

hclust_out_tr = cell(cdata.num_dsets,1);
fig_h_tr = figure;
sp_h_tr = cell(cdata.num_dsets,1);
hclust_out_cell = cell(cdata.num_dsets,1);
fig_h_cell = figure;
sp_h_cell = cell(cdata.num_dsets,1);
fig_ras = figure;
sp_ras = cell(cdata.num_dsets,1);
% fig_ras2 = figure;
% sp_ras2 = cell(cdata.num_dsets,1);

data_dim_est_full = cell(cdata.num_dsets,1);
ens_out_full = cell(cdata.num_dsets,1);

for n_dset = 1:cdata.num_dsets
    disp([dr_params.cond_name, ' dset ' num2str(n_dset)]);

    trial_types = cdata.trial_types_pr{n_dset};
    trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
    trial_data_sort_sm_pr = cdata.trial_data_sort_sm_pr{n_dset};
    
    %% select trials
    [tn_to_dred, trial_type_tag] = f_select_trial_type(dr_params.tt_to_dred_input, cdata, n_dset, ops);
    tt_to_dred = ops.context_types_all(tn_to_dred);
    
    trials_idx_dred = logical(sum(trial_types == tt_to_dred(:)' ,2));
    trial_peaks_dred = trial_peaks(:,trials_idx_dred);
    trial_types_dred = trial_types(trials_idx_dred);
    trial_data_sort_sm_pr = trial_data_sort_sm_pr(:,:,trials_idx_dred);
    
    %% dset specific params
    dr_params.n_dset = n_dset;
    dr_params.trial_type_tag = trial_type_tag;
    dr_params.tn_to_dred = tn_to_dred;
    dr_params.tt_to_dred = tt_to_dred;
    dr_params.volume_period = cdata.proc_data{n_dset}.frame_data.volume_period;
    dr_params.trial_t = cdata.trial_window_t{n_dset};
    dr_params.ctx_mmn = ops.context_types_all(cdata.ctx_mmn{n_dset});
    

    %% select responsive cells
    if ops.dred_params.use_responsive_cells
        resp_cells = cdata.peak_tuned_trials_combined{n_dset};
        resp_cells = logical(sum(resp_cells(:,tn_to_dred),2));
        trial_peaks_dred = trial_peaks_dred(resp_cells,:);
        trial_data_sort_sm_pr = trial_data_sort_sm_pr(resp_cells,:,:);
    end
    
    if sum(resp_cells)>5
        %%
        %f_hclust_estimate_num_clust(trial_peaks_dred, dr_params, ops)

        %% hclustering trials
        figure(fig_h_tr);
        sp_h_tr{n_dset} = subplot(3,5,n_dset);
        hc_params = dr_params;
        hc_params.subplot_ptr = sp_h_tr{n_dset};
        hc_params.method = ops.dred_params.hclust.method;
        hc_params.metric = ops.dred_params.hclust.plot_metric;
        hclust_out_tr{n_dset} = f_hcluster_trial2(trial_peaks_dred, trial_types_dred , hc_params, ops);
        dr_params.hclust_out_tr = hclust_out_tr{n_dset};
        
        %% hclustering cells 
        figure(fig_h_cell);
        sp_h_cell{n_dset} = subplot(3,5,n_dset);
        hc_params.subplot_ptr = sp_h_cell{n_dset};
        hclust_out_cell{n_dset} = f_hcluster_cell(trial_peaks_dred, trial_types_dred, hc_params, ops);
        dr_params.hclust_out_cell = hclust_out_cell{n_dset};
        %%
        %f_tsne(trial_peaks)


        %% compute dimensionality of full dsets
%         trial_peaks_dred_sort = trial_peaks_dred(:,hclust_out_tr{n_dset}.dend_order);
%         trial_peaks_dred_sort = trial_peaks_dred_sort(hclust_out_cell{n_dset}.dend_order,:);
%         trial_types_dred_sort = trial_types_dred(hclust_out_tr{n_dset}.dend_order);
        
        data_dim_est_full{n_dset} = f_ensemble_comp_data_dim2(trial_peaks_dred, 0);
        %% extract ensembles? 
        if ops.dred_params.do_ensamble_analysis
            params.cond_name = dr_params.cond_name;
            params.n_dset = dr_params.n_dset;
            params.normalize = 1;
            %params.num_comps = [];
            params.plot_stuff = 0;
            params.ensamble_method = ops.ensemb;

            ens_out_full{n_dset} = f_ensemble_analysis_peaks3(trial_peaks_dred, params, ops);
            
        end
        %%
        %ops.dred_params.hclust.sort_raster = 1;
        %raster_intput1 = trial_data_sort_sm_pr;
        raster_intput1 = trial_peaks_dred;
        dr_params.dend_order_cells = ens_out_full{n_dset}.cell_clust.dend_order;
        dr_params.clust_ident_cells = ens_out_full{n_dset}.cell_clust.clust_ident;
        dr_params.dend_order_trials = ens_out_full{n_dset}.trial_clust.dend_order;
        dr_params.clust_ident_trials = ens_out_full{n_dset}.trial_clust.clust_ident;
        figure(fig_ras);
        sp_ras{n_dset} = subplot(3,5,n_dset);
        f_hclust_raster(raster_intput1, trial_types_dred, sp_ras{n_dset}, dr_params, ops);
        
%         ops.dred_params.hclust.sort_raster = 0;
%         figure(fig_ras2);
%         sp_ras2{n_dset} = subplot(3,5,n_dset);
%         f_hclust_raster(trial_data_sort_sm_pr, trial_peaks_dred, trial_types_dred, sp_ras2{n_dset}, dr_params, ops);
%         
        
        
        %% distance metric
%         num_tt = numel(tn_to_dred);
%         if num_tt == 2
%             tt_index1 = trial_types_dred_sort == tn_to_dred(1);
%             f_ensemble_analysis_peaks(trial_peaks_dred_sort(:,tt_index1), trial_types_dred_sort(tt_index1), dr_params, ops);
%         end
        
        %% compute data dimensionality for cell range
        if ops.dred_params.do_dim_estimate
            dim_est_st = dr_params.dim_est_st;
            dd_idx = numel([dim_est_st.n_dset])+1;
            num_cells = size(trial_data_sort_sm_pr,1);
            interval1 = 5;
            num_repeats = 10;
            dd_cells_range = [interval1:interval1:num_cells num_cells];
            for n_cellr = 1:numel(dd_cells_range)
                for n_rep = 1:num_repeats
                    
                    samp_idx = randsample(num_cells, dd_cells_range(n_cellr));
                    data_dim_est = f_ensemble_comp_data_dim2(trial_peaks_dred_sort(samp_idx,:), 0);

                    %data_dim_est = f_ensemble_analysis_YS2(trial_data_sort_sm,trial_types_dred);
                    dim_est_st(dd_idx).cond_name = dr_params.cond_name;
                    dim_est_st(dd_idx).n_dset = n_dset;
                    dim_est_st(dd_idx).tt_to_dred = tt_to_dred;
                    dim_est_st(dd_idx).trial_type_tag = trial_type_tag;
                    dim_est_st(dd_idx).num_cells = num_cells;
                    dim_est_st(dd_idx).num_cells_samp = dd_cells_range(n_cellr);
                    dim_est_st(dd_idx).dimensionality_total = data_dim_est.dimensionality_total;
                    dim_est_st(dd_idx).dimensionality_corr = data_dim_est.dimensionality_corr;
                    dim_est_st(dd_idx).num_comp_est = data_dim_est.num_comps;
                    dim_est_st(dd_idx).d_explained = data_dim_est.d_explained;
                    dim_est_st(dd_idx).n_rep = n_rep;
                    dim_est_st(dd_idx).var_thresh_prc = data_dim_est.var_thresh_prc;             
                    dd_idx = dd_idx+1;
                end
            end
            dr_params.dim_est_st = dim_est_st;
        end
        
        %%
        

    end
    
end
figure(fig_h_tr);
suptitle(sprintf('%s; %s trials clust=%d; trials:%s', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, trial_type_tag));
figure(fig_h_cell);
suptitle(sprintf('%s; %s cells clust=%d; trials:%s', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, trial_type_tag));
figure(fig_ras);
suptitle(sprintf('%s; %s clust=%d; trials:%s sort', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, trial_type_tag));
% figure(fig_ras2);
% suptitle(sprintf('%s; %s clust=%d; trials:%s', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, trial_type_tag));

hdata_out.hclust_out_tr = hclust_out_tr;
hdata_out.hclust_out_cell = hclust_out_cell;
hdata_out.data_dim_est_full = data_dim_est_full;
end