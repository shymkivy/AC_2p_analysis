function f_dv_compute_stats2(app)

stat_window = [-2 4];
stat_resp_window = [.05 1];

%stat_resp_window = [-2 3];

z_thresh = app.ZthreshnewEditField.Value;

%%
n_pl = app.mplSpinner.Value;

num_cells = app.ddata.num_cells_pl{n_pl};
stim_times = app.ddata.stim_frame_index{n_pl};
%trig_window = app.working_ops.trial_num_baseline_resp_frames;
trial_types = app.ddata.trial_types{1};
MMN_freq = app.ddata.MMN_freq{1};
fr = 1000/double(app.ddata.proc_data{1}.frame_data.volume_period);

%%
stat_window_t = (ceil(stat_window(1)*fr):floor(stat_window(2)*fr))/fr;
stat_window_num_baseline_resp_frames = [sum(stat_window_t<=0) sum(stat_window_t>0)];   

%%
trial_data_sort = f_get_stim_trig_resp(app.cdata.S, stim_times, stat_window_num_baseline_resp_frames);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, MMN_freq, app.ops);

% get stim times for shuffle
if strcmpi(app.StatsourceDropDown.Value, 'All')
    pop_stim_times = stim_times;
    trial_data_sort_stat = trial_data_sort;
elseif strcmpi(app.StatsourceDropDown.Value, 'Freqs')
    stim_idx = logical(sum(trial_types == 1:app.ops.stim.num_freqs,2));
    pop_stim_times = stim_times(stim_idx);
    trial_data_sort_stat = trial_data_sort(:,:,stim_idx);
end

num_tt = numel(app.ops.context_types_all);
num_trials = numel(pop_stim_times);

pop_mean = cell(num_cells,1);
pop_z_factor = cell(num_cells,1);


wb = f_waitbar_initialize(app, 'Computing statistics...');
for n_cell = 1:num_cells
    f_waitbar_update(wb, n_cell/num_cells, 'Computing statistics...');
    
    resp_all = squeeze(trial_data_sort_stat(n_cell,:,:));
    
    if strcmpi(app.StatmethodDropDown.Value, 'Sample')
        
        samp_size = round(sum(sum(trial_types == 1:app.ops.stim.num_freqs,2))/app.ops.stim.num_freqs);
        num_samp = 1000;

        samp_mean = zeros(sum(stat_window_num_baseline_resp_frames), num_samp);
        for n_samp = 1:num_samp
            samp = randsample(num_trials, samp_size);
            samp_resp = resp_all(:,samp);
            samp_mean(:,n_samp) = mean(samp_resp,2);
        end

        resp_all_mean = mean(samp_mean,2);
        z_factor = std(samp_mean,[],2);
    elseif strcmpi(app.StatmethodDropDown.Value, 'Pop percentile')
        resp_all_mean = mean(resp_all,2);
        z_factor = (prctile(resp_all', 95)-resp_all_mean)/2;
    end
    pop_mean{n_cell} = resp_all_mean;
    pop_z_factor{n_cell} = z_factor;
end
f_waitbar_close(wb);

peak_in_win = false(num_cells, num_tt);
peak_is_sig = false(num_cells, num_tt);
cell_is_resp = false(num_cells, num_tt);
peak_val_all = zeros(num_cells, num_tt);
peak_t_all = zeros(num_cells, num_tt);
for n_cell = 1:num_cells
    for n_tt = 1:num_tt
        temp_resp = trial_data_sort_wctx(n_cell,:,trial_types_wctx == app.ops.context_types_all(n_tt));
        temp_resp2 = squeeze(temp_resp);
        
        mean_resp = mean(temp_resp2,2);
        
        [pks,locs] = findpeaks(mean_resp);
        [peak_val, max_idx2] = max(pks);
        max_idx = locs(max_idx2);
        
        peak_t_all(n_cell, n_tt) = stat_window_t(max_idx);
        peak_val_all(n_cell, n_tt) = peak_val;
        peak_in_win(n_cell, n_tt) = and(stat_window_t(max_idx)>=stat_resp_window(1),stat_window_t(max_idx)<=stat_resp_window(2));
        peak_is_sig(n_cell, n_tt) = peak_val>(pop_mean{n_cell}(max_idx)+pop_z_factor{n_cell}(max_idx)*z_thresh);
                
        if and(peak_in_win(n_cell, n_tt),peak_is_sig(n_cell, n_tt))
            cell_is_resp(n_cell, n_tt) = 1;
        end
        
    end
end

col1 = jet(10);

figure; hold on;
for n_tt = 1:10
    [f, x] = ksdensity(peak_t_all(peak_is_sig(:,n_tt),n_tt));
    plot(x, f, 'Color', col1(n_tt,:));
end

stats.pop_mean = pop_mean;
stats.pop_z_factor = pop_z_factor;
stats.cell_is_resp = cell_is_resp;
stats.peak_val_all = peak_val_all;
stats.peak_t_all = peak_t_all;
stats.stat_window_t = stat_window_t;
stats.z_thresh = z_thresh;

app.ddata.stats{n_pl} = stats;
app.data(app.current_data_idx,:).stats{n_pl} = stats;

end