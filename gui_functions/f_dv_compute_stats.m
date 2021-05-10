function f_dv_compute_stats(app)

stat_window = [-2 4];
stat_resp_window = [.05 1];

z_thresh = 2;

%%
n_pl = app.mplSpinner.Value;

num_cells = app.ddata.num_cells_pl{n_pl};
cuts_trace = logical(app.ddata.proc_data{1}.file_cuts_params{n_pl}.vid_cuts_trace);
num_t = numel(cuts_trace);
stim_times = app.ddata.stim_frame_index{n_pl};
%trig_window = app.working_ops.trial_num_baseline_resp_frames;
trial_types = app.ddata.trial_types{1};
fr = double(app.ddata.OA_data{n_pl}.ops.init_params_caiman.data.fr);

%%
stat_window_t = (ceil(stat_window(1)*fr):floor(stat_window(2)*fr))/fr;
stat_window_num_baseline_resp_frames = [sum(stat_window_t<=0) sum(stat_window_t>0)];   

%%
accepted_cells = app.ddata.OA_data{n_pl}.proc.comp_accepted;
cell_num_convert = find(accepted_cells);

%C = app.ddata.OA_data{n_pl}.est.C;
%Yra = app.ddata.OA_data{n_pl}.est.YrA;
if strcmpi(app.DeconvolutionmethodDropDown.Value, 'OA_deconv')
    S = app.ddata.OA_data{n_pl}.est.S;
elseif strcmpi(app.DeconvolutionmethodDropDown.Value, 'MCMC')
    %C = app.ddata.OA_data{n_pl}.proc.deconv.MCMC.C;
    S = app.ddata.OA_data{n_pl}.proc.deconv.MCMC.S;
elseif strcmpi(app.DeconvolutionmethodDropDown.Value, 'smooth_dfdt')
    S = app.ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.S;
end

if strcmpi(app.StatsourceDropDown.Value, 'All')
    pop_stim_times = stim_times;
elseif strcmpi(app.StatsourceDropDown.Value, 'Freqs')
    stim_idx = logical(sum(trial_types == 1:app.ops.stim.num_freqs,2));
    pop_stim_times = stim_times(stim_idx);
end


num_trials = numel(app.ops.context_types_all);

pop_mean = cell(num_cells,1);
pop_z_factor = cell(num_cells,1);
cell_is_resp = false(num_cells, num_trials);
peak_val_all = zeros(num_cells, num_trials);
peak_t_all = zeros(num_cells, num_trials);

wb = f_waitbar_initialize(app, 'Computing statistics...');
for n_cell = 1:num_cells
    f_waitbar_update(wb, n_cell/num_cells, 'Computing statistics...');
    %current_cell_raw = zeros(1,num_t);
    %current_cell_raw(cuts_trace) = C(cell_num_convert(n_cell),:) + Yra(cell_num_convert(n_cell),:);
    
    %current_cell_C = zeros(1,num_t);
    current_cell_spikes = zeros(1,num_t);
    
    if app.SmoothCheckBox.Value 
        S2 = f_smooth_gauss2(S(cell_num_convert(n_cell),:), app.SmoothsigmamsEditField.Value/1000*fr, 0);
    else
        S2 = S(cell_num_convert(n_cell),:);
    end
    
    %current_cell_C(cuts_trace) = C(cell_num_convert(n_cell),:);
    current_cell_spikes(cuts_trace) = S2;

    resp_all = squeeze(f_get_stim_trig_resp(current_cell_spikes, pop_stim_times, stat_window_num_baseline_resp_frames));

    if strcmpi(app.StatmethodDropDown.Value, 'Sample')
        num_trials = numel(pop_stim_times);
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
    
    for n_tr = 1:10%num_trials
        temp_resp = f_get_stim_trig_resp(current_cell_spikes, stim_times(trial_types == n_tr), stat_window_num_baseline_resp_frames);
        temp_resp2 = squeeze(temp_resp);
        
        mean_resp = mean(temp_resp2,2);
        
        [pks,locs] = findpeaks(mean_resp);
        [peak_val, max_idx2] = max(pks);
        max_idx = locs(max_idx2);
        
        peak_in_resp_win = and(stat_window_t(max_idx)>stat_resp_window(1),stat_window_t(max_idx)<stat_resp_window(2));
        peak_sig = peak_val>(resp_all_mean(max_idx)+z_factor(max_idx)*z_thresh);
        
        if and(peak_in_resp_win,peak_sig)
            cell_is_resp(n_cell, n_tr) = 1;
            peak_t_all(n_cell, n_tr) = stat_window_t(max_idx);
            peak_val_all(n_cell, n_tr) = peak_val;
        end
        
    end
end
f_waitbar_close(wb);

stats.pop_mean = pop_mean;
stats.pop_z_factor = pop_z_factor;
stats.cell_is_resp = cell_is_resp;
stats.peak_val_all = peak_val_all;
stats.peak_t_all = peak_t_all;
stats.stat_window_t = stat_window_t;

data_idx = strcmpi(app.data.experiment, app.ddata.experiment);
app.data.stats(data_idx) = stats;
app.ddata.stats = stats;
end