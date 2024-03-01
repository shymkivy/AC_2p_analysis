function f_dv_plot_pop_mmn(app)
plot_first_tr = 1;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = numel(data.dset_name);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);
num_t = sum(trial_frames);

if sum(strcmpi(app.trialtypeDropDown.Value, {'Context', 'Context_flip', 'Context_both_comb'}))
    tt_input = app.trialtypeDropDown.Value;
else
    tt_input = 'Context_both_comb';
end
tn_all = f_dv_get_trial_number(app, tt_input);
[num_gr, num_tn] = size(tn_all);

num_red = 10;
vol_per = app.ddata.proc_data{1}.frame_data.volume_period_ave;
plot_frames = round(num_red/vol_per*1000);
plot_t2 = ((1:plot_frames))*vol_per/1000;

params = f_dv_gather_params(app);

[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val2(app);
num_reg = size(region_num,1);

mouse_id = cell(num_gr, num_dsets);
dset_id = zeros(num_gr, num_dsets);

num_cells_all = zeros(num_dsets, num_gr, num_reg);
data_first_trial = cell(num_dsets, num_gr, num_reg);
data_rem_trials_mean = cell(num_dsets, num_gr, num_reg);
data_rem_trials_cellmean = cell(num_dsets, num_gr, num_reg);
data_rem_trials_all = cell(num_dsets, num_gr, num_reg);
data_mmn_all = cell(num_dsets, num_gr, num_reg);
fprintf('dset #/%d: ', num_dsets);
for n_dset = 1:num_dsets
    fprintf('..%d', n_dset);
    data1 = data(n_dset,:);
    stats1 = data1.stats{n_pl};
    params.n_dset = find(data1.idx == app.data.idx);
    cdata = f_dv_compute_cdata(data1, params);
    
    if app.ConverttoZCheckBox.Value
        st_mean_mean = stats1.stat_trials_mean_mean;
        st_mean_sem = stats1.stat_trials_mean_sem;
    else
        st_mean_mean = zeros(cdata.num_cells,1);
        st_mean_sem = ones(cdata.num_cells,1);
    end

    firing_rate = cdata.S_sm;
    mmn_freq = data1.MMN_freq{1};
    trial_types = data1.trial_types{1};
    trial_types_ctx = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);
    stim_times = data1.stim_frame_index{n_pl};
    num_tt_all = numel(trial_types);

    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
    trial_data_sort = (trial_data_sort - st_mean_mean)./st_mean_sem;

    reg_cell_labels = f_dv_get_area_label(app, data1);

    for n_gr = 1:num_gr
        mouse_id{n_gr, n_dset} = data1.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;

        tn1 = tn_all(n_gr, :);
        num_tn = numel(tn1);
        red_start = app.ops.context_types_all(tn1(1) - 7);
        red_end = app.ops.context_types_all(50 - tn1(3));
        mmn_start = find(trial_types == red_start);
        num_tr = numel(mmn_start);
        mmn_end = zeros(num_tr,1);

        for n_st = 1:num_tr
            curr_ntt = mmn_start(n_st);
            tt1 = trial_types(curr_ntt);
            while and(tt1 >= red_start, tt1 < red_end)
                curr_ntt = curr_ntt + 1;
                if curr_ntt > num_tt_all
                    tt1 = red_end;
                else
                    tt1 = trial_types(curr_ntt);
                end
            end
            mmn_end(n_st) = curr_ntt-1;
        end

        frame_start = stim_times(mmn_start);
        frame_end = stim_times(mmn_end) + trial_frames(2);
        stim_dur = frame_end - frame_start+1;
        max_stim_dur = max([max(stim_dur), plot_frames]);

        trial_data_sort_red = nan(cdata.num_cells, max_stim_dur, num_tr);
        
        for n_tr = 1:num_tr
            trial_data_sort_red(:, 1:stim_dur(n_tr), n_tr) = firing_rate(:,frame_start(n_tr):frame_end(n_tr));
        end
        
        trial_data_sort_red = (trial_data_sort_red - st_mean_mean)./st_mean_sem;
        
        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn1);
        resp_cells2 = logical(sum(resp_cells,2));
        
        for n_reg = 1:num_reg
            n_reg2 = region_num(n_reg,:);
            
            % get resp cells
            reg_cell_idx = logical(sum(reg_cell_labels == n_reg2,2));
            resp_cell_idx = logical(resp_cells2.*reg_cell_idx);
            
            num_cells = sum(resp_cell_idx);
            num_cells_all(n_dset, n_gr, n_reg) = num_cells;
            
            if num_cells

                trial_data_sort_red2 = trial_data_sort_red(resp_cell_idx,:,:);
                trial_data_sort2 = trial_data_sort(resp_cell_idx,:,:);
                
                data_first_trial{n_dset, n_gr, n_reg} = trial_data_sort_red2(:, 1:plot_frames, 1);
                data_rem_trials_mean{n_dset, n_gr, n_reg} = mean(trial_data_sort_red2(:, 1:plot_frames, 2:end), 3, "omitnan");
                data_rem_trials_cellmean{n_dset, n_gr, n_reg} = mean(permute(trial_data_sort_red2(:, 1:plot_frames, 2:end), [3 2 1]), 3, "omitnan");
                data_rem_trials_all{n_dset, n_gr, n_reg} = trial_data_sort_red2(:, 1:plot_frames, 2:end);
                
                data_mmn_all{n_dset, n_gr, n_reg} = zeros(num_cells, num_t, num_tn);
                for n_tn = 1:num_tn
                    ctx1 = app.ops.context_types_all(tn1(n_tn));

                    idx1 = logical(sum([trial_types, trial_types_ctx] == ctx1,2));
                    temp_resp = trial_data_sort2(:, :, idx1);
                    data_mmn_all{n_dset, n_gr, n_reg}(:,:,n_tn) = mean(temp_resp,3);
                end
            end
        end
    end
end
fprintf('\n');

num_effdsets = num_dsets*num_gr;

data_first_trial2 = reshape(data_first_trial, num_effdsets, num_reg);
data_rem_trials_mean2 = reshape(data_rem_trials_mean, num_effdsets, num_reg);
data_rem_trials_cellmean2 = reshape(data_rem_trials_cellmean, num_effdsets, num_reg);
data_rem_trials_all2 = reshape(data_rem_trials_all, num_effdsets, num_reg);
data_mmn_all2 = reshape(data_mmn_all, num_effdsets, num_reg);

mouse_id2 = reshape(mouse_id, num_effdsets, 1);
dset_id2 = reshape(dset_id, num_effdsets, 1);

[gr_id, gr_tag] = f_dv_combine_data(app, mouse_id2, dset_id2);

gr_all = unique(gr_id);
num_gr2 = numel(gr_all);

%%

transp = app.StimtranspEditField.Value;
freq_col = app.stimcolorSpinner.Value;

title_tag2 = sprintf('%s; %s; %s; %s reg', title_tag, app.ResponsivecellsselectDropDown.Value, tt_input, reg_tag);

n_gr = 1;
tn1 = tn_all(n_gr, :);
num_tn = numel(tn1);

for n_reg = 1:num_reg
    title_tag3 = sprintf('%s; %s', title_tag2, leg_list{n_reg});
    
    first_tr_all = cat(1, data_first_trial2{:, n_reg});
    tr_mean_all = cat(1, data_rem_trials_mean2{:, n_reg});

    data_mmn_all3 = cat(1, data_mmn_all2{:, n_reg});
    
    num_cells = size(first_tr_all,1);
    if num_cells
        
        first_tr_mean  = mean(first_tr_all,1, "omitnan");
        first_tr_sem = std(first_tr_all, "omitnan")/sqrt(num_cells-1);
        
        tr_mean_mean = mean(tr_mean_all, 1, "omitnan");
        tr_mean_sem = std(tr_mean_all, "omitnan")/sqrt(num_cells-1);

        data_mmn_all3_mean = mean(data_mmn_all3, 1);
        data_mmn_all3_sem = std(data_mmn_all3, [], 1)/sqrt(num_cells-1);
        
        if plot_first_tr
            ymin = min([first_tr_mean-first_tr_sem, tr_mean_mean-tr_mean_sem, data_mmn_all3_mean(:)'-data_mmn_all3_sem(:)'])*1.05;
            ymax = max([max([first_tr_mean+first_tr_sem, tr_mean_mean+tr_mean_sem, data_mmn_all3_mean(:)'+data_mmn_all3_sem(:)'])*1.05, app.maxYlimEditField.Value]);
        else
            ymin = min([tr_mean_mean-tr_mean_sem, data_mmn_all3_mean(:)'-data_mmn_all3_sem(:)'])*1.05;
            ymax = max([max([tr_mean_mean+tr_mean_sem, data_mmn_all3_mean(:)'+data_mmn_all3_sem(:)'])*1.05, app.maxYlimEditField.Value]);
        end

        color2 = app.ops.context_types_all_colors2{19};
        if plot_first_tr
            figure; axis tight; hold on
            if app.PlotstimCheckBox.Value
                for n_st = 1:num_red
                    r1 = rectangle('Position', [n_st-1 ymin 0.5 ymax-ymin]);
                    r1.FaceColor = [app.ops.context_types_all_colors2{freq_col} transp];
                    r1.EdgeColor = [app.ops.context_types_all_colors2{freq_col} transp];
                end
            end
            shadedErrorBar_YS(plot_t2, first_tr_mean, first_tr_sem, color2);
            ylim([ymin, ymax]);
            xlim([plot_t2(1), plot_t2(end)]);
            if app.ConverttoZCheckBox.Value
                ylabel('Z-score');
            else
                ylabel('response mag');
            end
            title(sprintf('first trial %s; %d cells', title_tag3, num_cells), 'interpreter', 'none')
        end
        figure; hold on;
        if app.PlotstimCheckBox.Value
            for n_st = 1:num_red
                r1 = rectangle('Position', [n_st-1 ymin 0.5 ymax-ymin]);
                r1.FaceColor = [app.ops.context_types_all_colors2{freq_col} transp];
                r1.EdgeColor = [app.ops.context_types_all_colors2{freq_col} transp];
            end
        end
        shadedErrorBar_YS(plot_t2, tr_mean_mean, tr_mean_sem, color2);
        ylim([ymin, ymax]);
        xlim([plot_t2(1), plot_t2(end)])
        if app.ConverttoZCheckBox.Value
            ylabel('Z-score');
        else
            ylabel('response mag');
        end
        title(sprintf('2-end trials %s; %d cells', title_tag3, num_cells), 'interpreter', 'none')

        figure; 
        for n_tn = 1:num_tn
            tt1 = tn1(n_tn);
            subplot(1, 3, n_tn); hold on;
            if app.PlotstimCheckBox.Value
                r1 = rectangle('Position', [0 ymin 0.5 ymax-ymin]);
                r1.FaceColor = [app.ops.context_types_all_colors2{freq_col} transp];
                r1.EdgeColor = [app.ops.context_types_all_colors2{freq_col} transp];
            end
            color2 = app.ops.context_types_all_colors2{tt1};
            shadedErrorBar_YS(plot_t, data_mmn_all3_mean(1,:,n_tn), data_mmn_all3_sem(1,:,n_tn), color2);
            ylim([ymin, ymax]);
            if app.ConverttoZCheckBox.Value
                ylabel('Z-score');
            else
                ylabel('response mag');
            end
        end
        sgtitle(sprintf('context trials %s; %d cells', title_tag3, num_cells), 'interpreter', 'none')
    end
end

%%

for n_reg = 1:num_reg
    title_tag3 = sprintf('%s; %s', title_tag2, leg_list{n_reg});

    data_rem_trials_cellmean3 = cat(1, data_rem_trials_cellmean2{:,n_reg});
    figure; 
    imagesc(data_rem_trials_cellmean3)
    title(sprintf('mean over cells; %s', title_tag3),  'interpreter', 'none')

    data_rem_trials_mean3 = cat(1, data_rem_trials_mean2{:,n_reg});
    figure; 
    imagesc(data_rem_trials_mean3)
    title(sprintf('mean over trials; %s', title_tag3),  'interpreter', 'none')
end


end