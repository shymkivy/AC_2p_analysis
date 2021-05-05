function f_dv_plot_freq_resp(app)

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
plot_t = app.ddata.trial_window{1}.trial_window_t;

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
    plot(plot_t, resp_all{n_fr}, 'color', [.6 .6 .6])
    plot(plot_t, mean(resp_all{n_fr},2), 'color', [1 0 1], 'LineWidth', 2)
    if rem(n_fr,5) ~= 1
        set(gca,'ytick',[])
    end
    title(sprintf('Freq %d', n_fr))
end
sgtitle(sprintf('Dset %s; Cell %d', app.ddata.experiment{1}, n_cell), 'Interpreter', 'none')
end