function f_dv_ensless_single_trial(app)

tn_all = f_dv_get_trial_number(app);

select_resp_cells = app.selectrespcellsCheckBox.Value;
resort_by_ens = app.resortbyensCheckBox.Value;
sort_trials = app.sorttrialsCheckBox.Value;
sort_with_full_firing_rate = app.sortwithfullfrCheckBox.Value;

ddata = app.ddata;
cdata = f_dv_get_cdata(app);
stats1 = ddata.stats{1};

num_cells = sum([cdata.num_cells]);
firing_rate = cat(1,cdata.S_sm);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    resp_cell = logical(sum(app.ddata.stats{n_pl}.resp_cells_peak(:,tn_all),2));
else
    resp_cells_all = cell(numel(app.ddata.stats),1);
    for n_pl2 = 1:numel(app.ddata.stats)
        resp_cells_all{n_pl2} = logical(sum(app.ddata.stats{n_pl2}.resp_cells_peak(:,tn_all),2));
    end
    resp_cell = cat(1,resp_cells_all{:});
end

resp_ens = logical(sum(app.ddata.ensemble_tuning_stats{1}.resp_cells_peak(:,tn_all),2));

ens_list = app.ddata.ensembles{1}.cells.ens_list(app.ddata.ensemble_stats{1}.accepted_ensembles);
resp_ens_list = ens_list(resp_ens);
ens_cells = false(num_cells, numel(resp_ens_list));
for n_ens = 1:numel(resp_ens_list)
    ens_cells(resp_ens_list{n_ens}, n_ens) = 1;
end

resp_all = [resp_cell, ens_cells];


%%
trial_types = app.ddata.trial_types{1};
stim_frame_index = app.ddata.stim_frame_index{1};

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[~, trial_frames] = f_dv_compute_window_t(trial_window, cdata.volume_period);

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_frames);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, app.ddata.MMN_freq{1}, app.ops);

tr_idx = logical(sum(trial_types_wctx == app.ops.context_types_all(tn_all)',2));
tr_data =trial_data_sort_wctx(:,:,tr_idx);

if select_resp_cells
    sel_resp_cells = logical(sum(resp_all,2));
    tr_data2 = tr_data(sel_resp_cells,:,:);
    resp_all2 = resp_all(sel_resp_cells,:);
    firing_rate2 = firing_rate(sel_resp_cells,:);
else
    tr_data2 = tr_data;
    resp_all2 = resp_all;
    firing_rate2 = firing_rate;
end

[num_cells2, ~, num_tr] = size(tr_data2);

hc_params.plot_dist_mat = 1;
hc_params.plot_clusters = 0;

if sort_trials
    tr_data_2d_tr = reshape(tr_data2, [], num_tr);
    hclust_out_trial = f_hcluster_wrap(tr_data_2d_tr', hc_params);
    
    tr_data3 = tr_data2(:,:,hclust_out_trial.dend_order);
else
    tr_data3 = tr_data2;
end

if sort_with_full_firing_rate
    hclust_out_cell = f_hcluster_wrap(firing_rate2, hc_params);
else
    tr_data_2d = reshape(tr_data3, num_cells2, []);
    hclust_out_cell = f_hcluster_wrap(tr_data_2d, hc_params);
end

tr_data_2d_sort = reshape(tr_data3(hclust_out_cell.dend_order,:,:), num_cells2, []);
resp_all_sort = resp_all2(hclust_out_cell.dend_order,:);

sort_ord = (1:num_cells2)';
if resort_by_ens
    temp_resp = resp_all_sort;
    for n_ens = 1:numel(resp_ens_list)
        [~, idx1] = sort(temp_resp(:,numel(resp_ens_list)-n_ens+2), 'descend');
        temp_resp = temp_resp(idx1,:);
        sort_ord = sort_ord(idx1);
    end
end

%%
figure; 
subplot(1,10,1:8);
imagesc(tr_data_2d_sort(sort_ord,:))
subplot(1,10, 9:10);
imagesc(resp_all_sort(sort_ord,:))

end