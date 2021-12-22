function f_dv_initialize(app)

app.regdatapathEditField.Value = app.gui_ops.reg_data_path;
app.matdatapathEditField.Value = app.gui_ops.mat_data_path;

app.gui_plots.A_image = imagesc(app.UIAxes, 0);
app.gui_plots.A_image.ButtonDownFcn = @(~,~) f_dv_button_down(app, app.gui_plots.A_image);
app.gui_plots.plot_current_contour = [];
axis(app.UIAxes, 'equal');
axis(app.UIAxes, 'tight');

hold(app.UIAxes2, 'on');
app.gui_plots.plot_raw = plot(app.UIAxes2, 0,0, 'color', [0.8500, 0.3250, 0.0980]);
app.gui_plots.plot_C = plot(app.UIAxes2, 0,0, 'color', [0, 0.4470, 0.7410]);
app.gui_plots.plot_spikes = plot(app.UIAxes2, 0,0, 'color', [0.4940, 0.1840, 0.5560]);
app.gui_plots.plot_stim_times = plot(app.UIAxes2, 0,0, 'color', [0, 0, 0]);
axis(app.UIAxes2, 'tight');
%axis(app.UIAxes2, 'tight');

app.gui_plots.freq_resp_fig = [];
app.gui_plots.ctx_resp_fig = [];
app.gui_plots.select_resp_fig = [];
app.gui_plots.contours_gobj = [];
app.gui_plots.registration_fig = [];

%% load default stats

app.newStatmethodDropDown.Items = {'shuff_pool', 'shuff_locwise', 'z_thresh'}; % unused
app.newStatmethodDropDown.Value = app.gui_ops.stats.stat_method;

app.newStatsourceDropDown.Items = {'Freqs_dd', 'Freqs', 'All'};
app.newStatsourceDropDown.Value = app.gui_ops.stats.stat_source;

app.newZthreshEditField.Value = app.gui_ops.stats.z_thresh;
app.newPvalEditField.Value = (1 - normcdf(app.newZthreshEditField.Value));
app.newPeakbintimesecEditField.Value = app.gui_ops.stats.peak_bin_time;
app.newNumshuffsampEditField.Value = app.gui_ops.stats.num_shuff_samp;
app.newbasewinEditField.Value = app.gui_ops.stats.base_resp_win(1);
app.newrespwinEditField.Value = app.gui_ops.stats.base_resp_win(2);
app.newLocothreshEditField.Value = app.gui_ops.stats.loco_thresh;

%%
app.plotfeatureDropDown.Items = {'peak loc', 'resp mag'};
app.plottypeDropDown.Items = {'kde', 'ecdf', 'histogram'};

app.gui_ops.save_var_list_pl = {'stats'};
app.gui_ops.save_var_list = {'data_dim_pca', 'data_dim_cv',...
                    'ensembles', 'ensemble_stats', 'ensemble_tuning_stats',...
                    'ensless_dim_est'};
                
                
app.RunallDropDown.Items = [app.gui_ops.save_var_list_pl app.gui_ops.save_var_list];

app.regiontoplotDropDown.Items = {'All', 'A1', 'A2', 'AAF', 'UF'};

app.DeconvolutionmethodDropDown.Items = {'smooth_dfdt'};

end