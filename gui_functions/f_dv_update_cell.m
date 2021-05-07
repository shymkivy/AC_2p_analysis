function f_dv_update_cell(app)

n_pl = app.mplSpinner.Value;

num_cells = app.ddata.num_cells_pl{n_pl};
app.CellSpinner.Value = min([app.CellSpinner.Value, num_cells]);
n_cell = app.CellSpinner.Value;

time1 = app.ddata.proc_data{1}.frame_data.frame_times_mpl{n_pl}/1000;

cuts_trace = logical(app.ddata.proc_data{1}.file_cuts_params{n_pl}.vid_cuts_trace);
num_t = numel(cuts_trace);

cell_num_conv = find(app.ddata.OA_data{n_pl}.proc.comp_accepted);

C = app.ddata.OA_data{n_pl}.est.C(cell_num_conv(n_cell),:);
Yra = app.ddata.OA_data{n_pl}.est.YrA(cell_num_conv(n_cell),:);
app.current_cell_raw = zeros(1,num_t);
app.current_cell_raw(cuts_trace) = C + Yra;

app.current_cell_C = zeros(1,num_t);
app.current_cell_spikes = zeros(1,num_t);

if strcmpi(app.DeconvolutionmethodDropDown.Value, 'OA_deconv')
    S = app.ddata.OA_data{n_pl}.est.S(cell_num_conv(n_cell),:);
elseif strcmpi(app.DeconvolutionmethodDropDown.Value, 'MCMC')
    C = app.ddata.OA_data{n_pl}.proc.deconv.MCMC.C{cell_num_conv(n_cell)};
    S = app.ddata.OA_data{n_pl}.proc.deconv.MCMC.S{cell_num_conv(n_cell)};
elseif strcmpi(app.DeconvolutionmethodDropDown.Value, 'smooth_dfdt')
    S = app.ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.S(cell_num_conv(n_cell),:);
end

if app.SmoothCheckBox.Value
    fr = double(app.ddata.OA_data{n_pl}.ops.init_params_caiman.data.fr);
    S = f_smooth_gauss2(S, app.SmoothsigmamsEditField.Value/1000*fr, 0);
end

app.current_cell_C(cuts_trace) = C;
app.current_cell_spikes(cuts_trace) = S;

if app.RawButton.Value
    app.gui_plots.plot_raw.XData = time1;
    app.gui_plots.plot_raw.YData = app.current_cell_raw;
else
    app.gui_plots.plot_raw.XData = 0;
    app.gui_plots.plot_raw.YData = 0;
end
if app.CButton.Value
    app.gui_plots.plot_C.XData = time1;
    app.gui_plots.plot_C.YData = app.current_cell_C;
else
    app.gui_plots.plot_C.XData = 0;
    app.gui_plots.plot_C.YData = 0;
end
if app.SpikesButton.Value
    app.gui_plots.plot_spikes.XData = time1;
    app.gui_plots.plot_spikes.YData = app.current_cell_spikes*app.ScaleEditField_sp.Value + app.ShiftEditField_sp.Value;
else
    app.gui_plots.plot_spikes.XData = 0;
    app.gui_plots.plot_spikes.YData = 0;
end

if app.StimtimesButton.Value
    if strcmpi(app.trialtypeDropDown.Value, 'all')
        idx_stim = app.ddata.stim_frame_index{n_pl};
    else
        idx_ctx = strcmpi(app.trialtypeDropDown.Value, app.ops.context_types_labels);
        tt = app.ops.context_types_all(idx_ctx);
        idx_stim = app.ddata.stim_frame_index{n_pl}(app.ddata.trial_types{1} == tt);
    end
    
    stim_trace = zeros(numel(time1),1);
    stim_trace(idx_stim) = 1;
    
    app.gui_plots.plot_stim_times.XData = time1;
    app.gui_plots.plot_stim_times.YData = stim_trace*app.ScaleEditField_st.Value+app.ShiftEditField_st.Value;
else
    app.gui_plots.plot_stim_times.XData = 0;
    app.gui_plots.plot_stim_times.YData = 0;
end

contours_accepted = app.ddata.OA_data{n_pl}.est.contours(app.ddata.OA_data{n_pl}.proc.comp_accepted);
temp_contours = contours_accepted{n_cell};

if isgraphics(app.gui_plots.plot_current_contour)
    delete(app.gui_plots.plot_current_contour);
end

hold(app.UIAxes, 'on');
app.gui_plots.plot_current_contour = plot(app.UIAxes, temp_contours(:,1), temp_contours(:,2), 'color', [0.75, 0, 0.75], 'LineWidth', 2);
hold(app.UIAxes, 'off');

SNR_accepted = app.ddata.OA_data{n_pl}.proc.SNR2_vals(app.ddata.OA_data{n_pl}.proc.comp_accepted);
app.SNREditField.Value = SNR_accepted(n_cell);

%%
firing_rate = app.current_cell_spikes;
stim_times = app.ddata.stim_frame_index{n_pl};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
trial_types = app.ddata.trial_types{1};

if strcmpi(app.StatsourceDropDown.Value, 'All')
    pop_stim_times = stim_times;
elseif strcmpi(app.StatsourceDropDown.Value, 'Freqs')
    stim_idx = logical(sum(trial_types == 1:app.ops.stim.num_freqs,2));
    pop_stim_times = stim_times(stim_idx);
end
app.stats.resp_all = squeeze(f_get_stim_trig_resp(firing_rate, pop_stim_times, trig_window));

if strcmpi(app.StatmethodDropDown.Value, 'Sample')
    num_trials = numel(pop_stim_times);
    samp_size = round(sum(sum(trial_types == 1:app.ops.stim.num_freqs,2))/app.ops.stim.num_freqs);
    num_samp = 1000;

    samp_mean = zeros(sum(trig_window), num_samp);
    for n_samp = 1:num_samp
        samp = randsample(num_trials, samp_size);
        samp_resp = app.stats.resp_all(:,samp);
        samp_mean(:,n_samp) = mean(samp_resp,2);
    end

    app.stats.resp_all_mean = mean(samp_mean,2);
    app.stats.z_factor = std(samp_mean,[],2);
elseif strcmpi(app.StatmethodDropDown.Value, 'Pop percentile')
    app.stats.resp_all_mean = mean(app.stats.resp_all,2);
    app.stats.z_factor = (prctile(app.stats.resp_all', 95)-app.stats.resp_all_mean)/2;
end

if app.UpdatefigsCheckBox.Value
    if isgraphics(app.gui_plots.freq_resp_fig)
        f_dv_plot_freq_resp(app)
    end
    if isgraphics(app.gui_plots.ctx_resp_fig)
        f_dv_plot_ctx_resp(app)
    end
    if isgraphics(app.gui_plots.select_resp_fig)
        f_dv_plot_select_trial(app)
    end
end

end