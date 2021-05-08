function f_dv_compute_stats(app)



n_pl = app.mplSpinner.Value;

num_cells = app.ddata.num_cells_pl{n_pl};
cuts_trace = logical(app.ddata.proc_data{1}.file_cuts_params{n_pl}.vid_cuts_trace);
num_t = numel(cuts_trace);
stim_times = app.ddata.stim_frame_index{n_pl};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
trial_types = app.ddata.trial_types{1};
fr = double(app.ddata.OA_data{n_pl}.ops.init_params_caiman.data.fr);

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


for n_cell = 1:num_cells
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

    resp_all = squeeze(f_get_stim_trig_resp(current_cell_spikes, pop_stim_times, trig_window));

    if strcmpi(app.StatmethodDropDown.Value, 'Sample')
        num_trials = numel(pop_stim_times);
        samp_size = round(sum(sum(trial_types == 1:app.ops.stim.num_freqs,2))/app.ops.stim.num_freqs);
        num_samp = 1000;

        samp_mean = zeros(sum(trig_window), num_samp);
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
    1
    %figure; plot()
end

end