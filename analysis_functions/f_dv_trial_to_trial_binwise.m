function f_dv_trial_to_trial_binwise(ddata, cdata0, params, ops)
% trial to trial analysis full firing rates

cdata = cat(1,cdata0{params.planes});
tn_all0 = f_dv_get_trial_number(params);
tn_all = tn_all0(:)';
%[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);
%[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);

stats1 = [ddata.stats{params.planes}];

resp_cells_all = cat(1, stats1.peak_resp_cells);
resp_cell = logical(sum(resp_cells_all(:,tn_all),2));

firing_rate = cat(1,cdata.S_sm);

% if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
%     resp_cell = logical(sum(app.ddata.stats{n_pl}.peak_resp_cells(:,tn_all),2));
% else
%     resp_cells_all = cell(numel(app.ddata.stats),1);
%     for n_pl2 = 1:numel(app.ddata.stats)
%         resp_cells_all{n_pl2} = logical(sum(app.ddata.stats{n_pl2}.peak_resp_cells(:,tn_all),2));
%     end
%     resp_cell = cat(1,resp_cells_all{:});
% end

%%

vol_per = mean(cat(1, cdata.volume_period));

trial_types = ddata.trial_types{1};
stim_frame_index = ddata.stim_frame_index{1};
[~, trial_frames] = f_dv_compute_window_t(params.trial_window, vol_per);

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_frames);
if ~isempty(ddata.MMN_freq{1})
    [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, ddata.MMN_freq{1}, ops);
else
    trial_data_sort_wctx = trial_data_sort;
    trial_types_wctx = trial_types;
end

tr_idx = logical(sum(trial_types_wctx == ops.context_types_all(tn_all)',2));
tr_loc_full = find(tr_idx);
tr_data = trial_data_sort_wctx(:,:,tr_idx);

resp_cells = f_dv_get_resp_vals_cells(stats1, tn_all(1,:), params);
resp_cells2 = logical(sum(resp_cells,2));

tr_data2 = tr_data(resp_cells2,:,:);
resp_all2 = resp_cell(resp_cells2,:);
firing_rate2 = firing_rate(resp_cells2,:);

[num_cells2, ~, num_tr] = size(tr_data2);

hc_params.plot_dist_mat = 1;
hc_params.plot_clusters = 0;
hc_params.num_clust = 1;

if params.sort_trials
    tr_data_2d_tr = reshape(tr_data2, [], num_tr);
    hclust_out_trial = f_hcluster_wrap(tr_data_2d_tr', hc_params);
    tr_data3 = tr_data2(:,:,hclust_out_trial.dend_order);
    
    SI = 1-hclust_out_trial.dist;
    SI_vals = tril(SI,-1);
    SI_vals(SI_vals==0) = [];

    figure; 
    im1 = imagesc(1-hclust_out_trial.dist); axis equal tight;
    title(sprintf('%s, trial %d, Mean corr = %.3f', ddata.dset_name_full{1}, tn_all, mean(SI_vals)), 'interpreter', 'none');
    im1.Parent.CLim(2) = 0.6;
    colorbar;
else
    tr_data3 = tr_data2;
end

% figure; plot(hclust_out_trial.clust_ident(hclust_out_trial.dend_order))
% x = tr_loc_full(hclust_out_trial.clust_ident==2)-1;
% x(x==0) = [];
% figure; histogram(trial_types(x));
% 

if params.sort_with_full_firing_rate
    hclust_out_cell = f_hcluster_wrap(firing_rate2, hc_params);
else
    tr_data_2d = reshape(tr_data3, num_cells2, []);
    hclust_out_cell = f_hcluster_wrap(tr_data_2d, hc_params);
end

tr_data_2d_sort = reshape(tr_data3(hclust_out_cell.dend_order,:,:), num_cells2, []);
resp_all_sort = resp_all2(hclust_out_cell.dend_order,:);

sort_ord = (1:num_cells2)';
if params.resort_by_ens
    temp_resp = resp_all_sort;
    for n_ens = 1:numel(resp_ens_list)
        [~, idx1] = sort(temp_resp(:,numel(resp_ens_list)-n_ens+2), 'descend');
        temp_resp = temp_resp(idx1,:);
        sort_ord = sort_ord(idx1);
    end
end

figure; 
subplot(1,10,1:8);
imagesc(tr_data_2d_sort(sort_ord,:))
subplot(1,10, 9:10);
imagesc(resp_all_sort(sort_ord,:))
end