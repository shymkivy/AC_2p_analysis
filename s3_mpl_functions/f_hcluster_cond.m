function [dr_params, hdata_out] = f_hcluster_cond(cdata, dr_params, ops)

num_dsets = numel(cdata.area);

hclust_out_tr = cell(num_dsets,1);
if ops.dred_params.hclust.plot_hclust_trials
    fig_h_tr = figure;
    sp_h_tr = cell(num_dsets,1);
end

hclust_out_cell = cell(num_dsets,1);
if ops.dred_params.hclust.plot_hclust_cells
    fig_h_cell = figure;
    sp_h_cell = cell(num_dsets,1);
end
fig_ras = figure;
sp_ras = cell(num_dsets,1);
% fig_ras2 = figure;
% sp_ras2 = cell(num_dsets,1);

data_dim_est_full = cell(num_dsets,1);
ens_out_full = cell(num_dsets,1);

for n_dset = 1:num_dsets
    disp([dr_params.cond_name, ' dset ' num2str(n_dset)]);

    trial_types = cdata.trial_types_wctx{n_dset};
    trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
    trial_data_sort_sm_wctx = cdata.trial_data_sort_sm_wctx{n_dset};
    
    %% select trials
    [tn_to_dred, trial_type_tag] = f_select_trial_type(dr_params.tt_to_dred_input, cdata, n_dset, ops);
    tt_to_dred = ops.context_types_all(tn_to_dred);
    
    trials_idx_dred = logical(sum(trial_types == tt_to_dred(:)' ,2));
    trial_peaks_dred = trial_peaks(:,trials_idx_dred);
    trial_types_dred = trial_types(trials_idx_dred);
    trial_data_sort_sm_wctx = trial_data_sort_sm_wctx(:,:,trials_idx_dred);
    
    %% dset specific params
    dr_params.n_dset = n_dset;
    dr_params.trial_type_tag = trial_type_tag;
    dr_params.tn_to_dred = tn_to_dred;
    dr_params.tt_to_dred = tt_to_dred;
    dr_params.volume_period = cdata.proc_data{n_dset}.frame_data.volume_period;
    dr_params.trial_t = cdata.trial_window{n_dset}.trial_window_t;
    dr_params.ctx_mmn = ops.context_types_all(cdata.ctx_mmn{n_dset});
    

    %% select responsive cells
    if ops.dred_params.use_responsive_cells
        resp_cells = cdata.peak_tuned_trials_combined{n_dset};
        resp_cells = logical(sum(resp_cells(:,tn_to_dred),2));
        trial_peaks_dred = trial_peaks_dred(resp_cells,:);
        trial_data_sort_sm_wctx = trial_data_sort_sm_wctx(resp_cells,:,:);
    end
    
    if sum(resp_cells)>5
        %%
        %f_hclust_estimate_num_clust(trial_peaks_dred, dr_params, ops)

        %% hclustering trials
        hc_params = dr_params;
        if ops.dred_params.hclust.plot_hclust_trials
            figure(fig_h_tr);
            sp_h_tr{n_dset} = subplot(3,5,n_dset);
            hc_params.subplot_ptr = sp_h_tr{n_dset};
        end
        hc_params.method = ops.dred_params.hclust.method;
        hc_params.metric = ops.dred_params.hclust.plot_metric;
        hc_params.plot_dist_mat = ops.dred_params.hclust.plot_hclust_trials;
        hc_params.plot_clusters = 0;
        hc_params_cells = hc_params;
        
        hc_params.sample_types = f_tt_to_tn(trial_types_dred, ops);
        hc_params.sample_types_colors = ops.context_types_all_colors2;
        %hc_params.XY_label = 'Trials';
        
        hclust_out_tr{n_dset} = f_hcluster_wrap(trial_peaks_dred', hc_params);
        dr_params.hclust_out_tr = hclust_out_tr{n_dset};
        %% hclustering cells 
        if ops.dred_params.hclust.plot_hclust_cells
            figure(fig_h_cell);
            sp_h_cell{n_dset} = subplot(3,5,n_dset);
            hc_params.subplot_ptr = sp_h_cell{n_dset};
        end
        hc_params_cells.plot_dist_mat = ops.dred_params.hclust.plot_hclust_cells;
        %hc_params_cells.XY_label = 'Cells';
        
        hclust_out_cell{n_dset} = f_hcluster_wrap(trial_peaks_dred, hc_params_cells);
        dr_params.hclust_out_cell = hclust_out_cell{n_dset};
        %%
        %f_tsne(trial_peaks)
        
        %%
%         trial_peaks_dred_sort = trial_peaks_dred(:,hclust_out_tr{n_dset}.dend_order);
%         trial_peaks_dred_sort = trial_peaks_dred_sort(hclust_out_cell{n_dset}.dend_order,:);
%         trial_types_dred_sort = trial_types_dred(hclust_out_tr{n_dset}.dend_order);

        %% compute dimensionality of full dsets

        params_dd.total_dim_thresh = ops.ensemb.total_dim_thresh;
        params_dd.corr_comp_thresh = ops.ensemb.corr_comp_thresh;
        params_dd.normalize = ops.ensemb.normalize; %; %'norm_full' 'norm_mean' 'none'
        data_dim_est_full{n_dset} = f_ensemble_comp_data_dim2(trial_peaks_dred, params_dd);
        %% extract ensembles? 
        if ops.dred_params.do_ensamble_analysis
            params_ens.cond_name = dr_params.cond_name;
            params_ens.n_dset = dr_params.n_dset;
            params_ens.ensamble_method = ops.ensemb.method;
            params_ens.corr_comp_thresh = ops.ensemb.corr_comp_thresh;
            params_ens.normalize = ops.ensemb.normalize; %; %'norm_full' 'norm_mean' 'none'
            params_ens.num_comps = [];
            params_ens.plot_stuff = 0;
            params_ens.use_LR_proj = 0;
            params_ens.ensamble_extraction = 'thresh'; % 'thresh', 'clust'
            params_ens.ensamble_extraction_thresh = 'signal_clust_thresh'; % 'shuff'. 'clust_thresh', 'signal_z'
            params_ens.hcluster_method = 'average';
            params_ens.hcluster_distance_metric = 'cosine';
            ens_out_full{n_dset} = f_ensemble_analysis_peaks3(trial_peaks_dred, params_ens);
        end
        %%
        %ops.dred_params.hclust.sort_raster = 1;
        raster_plot_intput1 = trial_data_sort_sm_wctx;
        %raster_plot_intput1 = trial_peaks_dred;
        trial_types_input1 = trial_types_dred;
        
        dr_params.clust_ident_cells = ens_out_full{n_dset}.cells.clust_ident;
        dr_params.clust_ident_trials = ens_out_full{n_dset}.trials.clust_ident;
        
        dr_params.dend_order_cells = ens_out_full{n_dset}.cells.dend_order;
        dr_params.dend_order_trials = ens_out_full{n_dset}.trials.dend_order;
        %dr_params.dend_order_cells = dr_params.hclust_out_cell.dend_order;
        %dr_params.dend_order_trials = dr_params.hclust_out_tr.dend_order;
        dr_params.dim_corr = ens_out_full{n_dset}.data_dim_est.dimensionality_corr;
        figure(fig_ras);
        sp_ras{n_dset} = subplot(3,5,n_dset);
        f_hclust_raster(raster_plot_intput1, trial_types_input1, sp_ras{n_dset}, dr_params, ops);
        
%         ops.dred_params.hclust.sort_raster = 0;
%         figure(fig_ras2);
%         sp_ras2{n_dset} = subplot(3,5,n_dset);
%         f_hclust_raster(trial_data_sort_sm_wctx, trial_peaks_dred, trial_types_dred, sp_ras2{n_dset}, dr_params, ops);
%         
        
        
        %% distance metric
%         num_tt = numel(tn_to_dred);
%         if num_tt == 2
%             f_dist_metric(trial_peaks_dred, trial_types_dred, dr_params)
%         end
        
        %% compute data dimensionality for cell range
        if ops.dred_params.do_dim_estimate
            dim_est_st = dr_params.dim_est_st;
            dd_idx = numel([dim_est_st.n_dset])+1;
            num_cells = size(trial_data_sort_sm_wctx,1);
            interval1 = 3;
            num_repeats = 10;
            dd_cells_range = interval1:interval1:min(num_cells,50);
            for n_cellr = 1:numel(dd_cells_range)
                for n_rep = 1:num_repeats
                    
                    samp_idx = randsample(num_cells, dd_cells_range(n_cellr));
                    data_dim_est = f_ensemble_comp_data_dim2(trial_peaks_dred(samp_idx,:), params_dd);
                    %data_dim_est = f_ensemble_analysis_YS2(trial_data_sort_sm,trial_types_dred);
                    dim_est_st(dd_idx).cond_name = dr_params.cond_name;
                    dim_est_st(dd_idx).n_dset = n_dset;
                    dim_est_st(dd_idx).tt_to_dred = tt_to_dred;
                    dim_est_st(dd_idx).trial_type_tag = trial_type_tag;
                    dim_est_st(dd_idx).num_cells = num_cells;
                    dim_est_st(dd_idx).num_cells_samp = dd_cells_range(n_cellr);
                    dim_est_st(dd_idx).dimensionality_total = data_dim_est.dimensionality_total;
                    dim_est_st(dd_idx).dimensionality_total_norm = data_dim_est.dimensionality_total_norm;
                    dim_est_st(dd_idx).dimensionality_total_norm_shuff = data_dim_est.dimensionality_total_norm_shuff;
                    dim_est_st(dd_idx).dimensionality_first_comp_size = data_dim_est.dimensionality_first_comp_size;
                    dim_est_st(dd_idx).dimensionality_corr = data_dim_est.dimensionality_corr;
                    dim_est_st(dd_idx).num_comp_est = data_dim_est.num_comps;
                    dim_est_st(dd_idx).d_explained = data_dim_est.d_explained;
                    dim_est_st(dd_idx).n_rep = n_rep;          
                    dd_idx = dd_idx+1;
                end
            end
            dr_params.dim_est_st = dim_est_st;
        end
        
        %%
        

    end
    
end
if ops.dred_params.hclust.plot_hclust_trials
    figure(fig_h_tr);
    suptitle(sprintf('%s; %s trials clust=%d; trials:%s', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, trial_type_tag));
end
if ops.dred_params.hclust.plot_hclust_cells
    figure(fig_h_cell);
    suptitle(sprintf('%s; %s cells clust=%d; trials:%s', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, trial_type_tag));
end
figure(fig_ras);
suptitle(sprintf('%s; %s clust=%d; trials:%s sort', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, trial_type_tag));
% figure(fig_ras2);
% suptitle(sprintf('%s; %s clust=%d; trials:%s', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, trial_type_tag));

hdata_out.hclust_out_tr = hclust_out_tr;
hdata_out.hclust_out_cell = hclust_out_cell;
hdata_out.data_dim_est_full = data_dim_est_full;
end