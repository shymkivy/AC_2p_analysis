function f_dv_plot_tuning(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;
ddata = app.ddata;

firing_rate = app.current_cell_spikes;
trial_types = app.ddata.trial_types{1};
stim_times = app.ddata.stim_frame_index{n_pl};

stats1 = app.ddata.stats{n_pl};

if app.ConverttoZCheckBox.Value
    st_mean_mean = stats1.stat_trials_mean_mean(n_cell);
    st_mean_sem = stats1.stat_trials_mean_sem(n_cell);
else
    st_mean_mean = 0;
    st_mean_sem = 1;
end

%grating_angles = data.proc_data{1,1}.stim_params.grating_angles;

tn_all = 1:10;

if strcmpi(app.TuningfreatureDropDown.Value, 'peaks')
    stats2 = stats1.peak;
elseif strcmpi(app.TuningfreatureDropDown.Value, 'onset')
    stats2 = stats1.onset;
elseif strcmpi(app.TuningfreatureDropDown.Value, 'offset')
    stats2 = stats1.offset;
end

figure;
polarplot(grating_angles*2 - pi/2, 1:10)

app.PhaseplotCheckBox.Value



end