function f_dv_plot_pop_ave_freq(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = numel(data.experiment);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);
num_t = sum(trial_frames);

tn_all = f_dv_get_trial_number2(app);
[num_tn_gr, num_tn] = size(tn_all);

params = f_dv_gather_params(app);

[region_num_all, reg_tag, leg_list] = f_dv_get_region_sel_val(app);
num_reg = size(region_num_all,1);

mouse_id = cell(num_tn_gr, num_dsets);
dset_id = zeros(num_tn_gr, num_dsets);

resp_all_mean = cell(num_tn_gr, num_dsets, num_reg);
resp_all_trials = cell(num_tn_gr, num_dsets, num_reg, num_tn);
resp_num_trials = zeros(num_tn_gr, num_dsets, num_reg, num_tn);
cell_counts = zeros(num_tn_gr, num_dsets, num_reg);

fprintf('dset #/%d: ', num_dsets);
for n_dset = 1:num_dsets
    fprintf('..%d', n_dset);
    data1 =  data(n_dset,:);
    stats1 = data1.stats{n_pl};
    params.n_dset = find(data1.idx == app.data.idx);

    cdata = f_dv_compute_cdata(data1, params);

    firing_rate = cdata.S_sm;
    mmn_freq = data1.MMN_freq{1};
    trial_types = data1.trial_types{1};
    trial_types_ctx = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);

    stim_times = data1.stim_frame_index{n_pl};
    
    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);

    reg_cell_labels = f_dv_get_area_label(app, data1);
    
    for n_tngr = 1:num_tn_gr
        mouse_id{n_tngr, n_dset} = data1.mouse_id{1};
        dset_id(n_tngr, n_dset) = n_dset;

        tn1 = f_dv_proc_tn(tn_all(n_tngr,:), mmn_freq);

        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn1);
        resp_cells2 = sum(resp_cells,2);
        %cell_is_resp = stats1.peak_resp_cells(:,tn_all);

        for n_reg = 1:num_reg
            region_num = region_num_all(n_reg,:);
            
            % get resp cells
            reg_cell_idx = logical(sum(reg_cell_labels == region_num,2));
            resp_cell_idx = logical(resp_cells2.*reg_cell_idx);

            num_cells = sum(resp_cell_idx);

            if num_cells
                cell_counts(n_tngr, n_dset, n_reg) = num_cells;
                resp_all_mean{n_tngr, n_dset, n_reg} = zeros(num_cells, num_t, num_tn);
                for n_tn = 1:num_tn
                    ctx1 = app.ops.context_types_all(tn1(n_tn));
                    idx1 = logical(sum([trial_types, trial_types_ctx] == ctx1,2));
                    temp_resp = trial_data_sort(resp_cell_idx,:,idx1);
                    resp_all_mean{n_tngr, n_dset, n_reg}(:,:,n_tn) = mean(temp_resp,3);
                    resp_all_trials{n_tngr, n_dset, n_reg, n_tn} = temp_resp;
                    resp_num_trials(n_tngr, n_dset, n_reg, n_tn) = size(temp_resp,3);
                end
            end
        end
    end
end
fprintf('\n');

num_effdsets = num_dsets*num_tn_gr;

resp_all_mean2 = reshape(resp_all_mean, num_effdsets, num_reg);
resp_all_trials2 = reshape(resp_all_trials, num_effdsets, num_reg, num_tn);
resp_num_trials2 = reshape(resp_num_trials, num_effdsets, num_reg, num_tn);
cell_counts2  = reshape(cell_counts, num_effdsets, num_reg);
mouse_id2 = reshape(mouse_id, num_effdsets, 1);
dset_id2 = reshape(dset_id, num_effdsets, 1);

[gr_id, gr_tag] = f_dv_combine_data(app, mouse_id2, dset_id2);


%%
reg_all = app.ops.regions_to_analyze;

ymin = 0;
ymax = 0;

title_tag2 = sprintf('%s; %s reg', title_tag, reg_tag);


for n_reg = 1:num_reg
    if cell_counts2(n_reg)
        resp_all2 = cat(1, resp_all_mean2{:, n_reg});
        tr_mean = squeeze(mean(resp_all2,1));
        ymin = min([min(tr_mean(:)), ymin]);
        ymax = max([max(tr_mean(:)), ymax]);
    end
end

for n_reg = 1:num_reg
    if cell_counts2(n_reg)
        num_tn2 = sum(tn_all < 11);
        if num_tn2
            figure();
            resp_all2 = cat(1, resp_all_mean2{:, n_reg});
            tr_mean = squeeze(mean(resp_all2,1));
            tr_std = squeeze(std(resp_all2, [], 1));
            
            for n_tn1 = 1:10
                n_tn = find(tn_all == n_tn1);
                subplot(1, 10, n_tn);
                plot(tr_mean(:, n_tn));
                ylim([ymin ymax]);
                title(sprintf('tn%d', tn_all2(n_tn1)));
            end
            sgtitle(sprintf('Freqs %s; %s', title_tag2, reg_all{n_reg}), 'interpreter', 'none')
        end
        
        tn2_idx = tn_all > 10;
        num_tn2 = sum(tn2_idx);
        tn_all2 = tn_all(tn2_idx);
        if num_tn2
            figure();
            resp_all2 = cat(1, resp_all_mean2{:, n_reg});
            tr_mean = squeeze(mean(resp_all2,1));
            tr_std = squeeze(std(resp_all2, [], 1));

            for n_tn1 = 1:num_tn2
                n_tn = find(tn_all == tn_all2(n_tn1));
                subplot(1, num_tn2, n_tn1);
                plot(tr_mean(:, n_tn));
                ylim([ymin ymax]);
                title(sprintf('tn%d', tn_all2(n_tn1)));
            end
            sgtitle(sprintf('Context %s; %s', title_tag2, reg_all{n_reg}), 'interpreter', 'none')
        end
    end
end



end