function f_dv_update_dset_info(app)

app.mplSpinner.Value = min([app.mplSpinner.Value, app.ddata.num_planes]);
app.MMNfreqEditField.Value = num2str(app.ddata.MMN_freq{1});

%%
app.FramerateEditField.Value = 1000/app.ddata.proc_data{1}.frame_data.volume_period;
frame_period = 1/app.FramerateEditField.Value;
trial_window = app.working_ops.trial_window;
trial_window_t = (ceil(trial_window(1)/frame_period):floor(trial_window(2)/frame_period))*frame_period;

app.working_ops.trial_window_t = trial_window_t;
app.working_ops.trial_num_baseline_resp_frames = [sum(trial_window_t<=0) sum(trial_window_t>0)];     

%%
n_pl = app.mplSpinner.Value;

deconv_methods = {};
if ~isempty(app.ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.S)
    deconv_methods = [deconv_methods, {'smooth_dfdt'}];
end
if ~isempty(cat(1,app.ddata.OA_data{n_pl}.proc.deconv.MCMC.S{:}))
    deconv_methods = [deconv_methods, {'MCMC'}];
end
if ~isempty(app.ddata.OA_data{n_pl}.est.S)
    deconv_methods = [deconv_methods, {'OA_deconv'}];
end
if ~isempty(cat(1,app.ddata.OA_data{n_pl}.proc.deconv.c_foopsi.S{:}))
    deconv_methods = [deconv_methods, {'constrained_foopsi'}];
end

app.DeconvolutionmethodDropDown.Items = deconv_methods;

f_dv_compute_stats(app);
f_dv_update_A(app);
f_dv_update_cell(app);


end