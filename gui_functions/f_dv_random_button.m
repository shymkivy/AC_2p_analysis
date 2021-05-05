function f_dv_random_button(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

if n_pl>1
   num_cells_pl2 = app.ddata.num_cells_pl(1:(n_pl-1));
   n_cell_full = sum(cat(1,num_cells_pl2{:})) + n_cell;
else
    n_cell_full = n_cell;
end

firing_rate = app.ddata.firing_rate_smooth{n_pl};
trial_types = app.ddata.trial_types_wctx{1};
stim_times = app.ddata.stim_frame_index{n_pl};
trig_window = app.ddata.trial_window{1}.trial_num_baseline_resp_frames;

resp_all = cell(10,1);
y_lim_max = 0;
y_lim_min = 0;
for n_fr = 1:10
    temp_resp = f_get_stim_trig_resp(firing_rate(n_cell,:), stim_times(trial_types == n_fr), trig_window);
    resp_all{n_fr} = squeeze(temp_resp);
    y_lim_max = max([y_lim_max max(resp_all{n_fr}(:))]);
    y_lim_min = min([y_lim_min min(resp_all{n_fr}(:))]);
end

figure; 
for n_fr = 1:10
    subplot(2,5,n_fr); hold on; axis tight; ylim([y_lim_min, y_lim_max]);
    plot(resp_all{n_fr}, 'color', [.6 .6 .6])
    plot(mean(resp_all{n_fr},2), 'color', [1 0 1], 'LineWidth', 2)
    title(sprintf('Freq %d', n_fr))
end

trial_window_t = app.ddata.trial_window{1}.trial_window_t;




resp = f_get_stim_trig_resp(firing_rate(n_cell,:), stim_times, trig_window);

resp1 = squeeze(resp);

trial_ave_z = app.ddata.trial_ave_z{1};
trial_ave_z_cell = trial_ave_z(n_cell_full,:,:);

if strcmpi(app.trialtypeDropDown.Value, 'all')
    
else
    ctx_idx = strcmpi(app.trialtypeDropDown.Value, app.ops.context_types_labels);
end

tt = app.ops.context_types_all(ctx_idx);

idx_trials = app.ddata.trial_types{1} == tt;
ops1 = app.ops;
ops1.plot_combined = 1;
f_mpl_plot_ctx_cell(resp1,trial_window_t,ops1)

figure; hold on;
plot(resp1, 'color', [.6 .6 .6])
plot(mean(resp1,2), 'color', [1 0 1], 'LineWidth', 2)
title(sprintf('%s', app.trialtypeDropDown.Value))


end