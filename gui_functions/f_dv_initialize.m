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

%%
gui_defaults = app.gui_ops.gui_defaults;
ops = app.gui_ops.ops;
stats = app.gui_ops.stats;
est_params_pca = app.gui_ops.est_params_pca;
est_params_cv = app.gui_ops.est_params_cv;
ens_params = app.gui_ops.ens_params;

%% gui defaults
%% general
app.plot_BaserespwinEditField.Value = f_array_to_str(ops.plot_window);
app.analysis_BaserespwinEditField.Value = f_array_to_str(ops.analysis_window);

%%
app.stats_StatmethodDropDown.Items = gui_defaults.stat_method_options;
app.stats_StatsourceDropDown.Items = gui_defaults.stat_source_options;

app.dimestpca_normalizeDropDown.Items = gui_defaults.normalize_options;
app.dimestpca_shufflemethodDropDown.Items = gui_defaults.shuffle_method_options;

app.dimestcv_normalizeDropDown.Items = gui_defaults.normalize_options;
app.dimestcv_ensmethodDropDown.Items = gui_defaults.est_dim_cv_method_options;

app.ens_ensemethodDropDown.Items = gui_defaults.ens_method_options;
app.ens_normalizeDropDown.Items = gui_defaults.normalize_options;
app.ens_ensextractionDropDown.Items = gui_defaults.ens_extraction_options;
app.ens_extractionthreshDropDown.Items = gui_defaults.ens_extraction_thresh_options;
app.ens_hclustmethodDropDown.Items = gui_defaults.hcluster_method_options;
app.ens_hclustmetricDropDown.Items = gui_defaults.hcluster_distance_metric_options;

%% load default stats

app.stats_StatmethodDropDown.Value = stats.stat_method;
app.stats_StatsourceDropDown.Value = stats.stat_source;
app.stats_ZthreshEditField.Value = stats.z_thresh;
app.stats_PvalEditField.Value = (1 - normcdf(app.stats_ZthreshEditField.Value));
app.stats_PeakbintimesecEditField.Value = stats.peak_bin_time;
app.stats_NumshuffsampEditField.Value = stats.num_shuff_samp;
app.stats_BaserespwinEditField.Value = f_array_to_str(stats.base_resp_win);
app.stats_LocothreshEditField.Value = stats.loco_thresh;

%% load default dim est pca 

app.dimestpca_normalizeDropDown.Value = est_params_pca.normalize;
app.dimestpca_shufflemethodDropDown.Value = est_params_pca.shuffle_method;
app.dimestpca_totaldimthreshEditField.Value = est_params_pca.total_dim_thresh;
app.dimestpca_dimestnumrepsEditField.Value = est_params_pca.dim_est_num_reps;
app.dimestpca_PlotstuffCheckBox.Value = est_params_pca.plot_stuff;

%% load default dim est cv

app.dimestcv_normalizeDropDown.Value = est_params_cv.normalize;
app.dimestcv_ensmethodDropDown.Value = est_params_cv.ensamble_method;
app.dimestcv_shuffledatachunksCheckBox.Value = est_params_cv.shuffle_data_chunks;

app.dimestcv_sm_centerEditField.Value = est_params_cv.smooth_SD_center;
app.dimestcv_sm_rangeEditField.Value = est_params_cv.smooth_SD_range;
app.dimestcv_sm_countEditField.Value = est_params_cv.smooth_SD_count;
app.dimestcv_smoothSDwindowEditField.Value = num2str(f_dv_get_dim_est_cv_sm(app));

app.dimestcv_nc_centerEditField.Value = est_params_cv.num_comp_center;
app.dimestcv_nc_rangeEditField.Value = est_params_cv.num_comp_range;
app.dimestcv_nc_countEditField.Value = est_params_cv.num_comp_count;
app.dimestcv_numcompwindowEditField.Value = num2str(f_dv_get_dim_est_cv_nc(app));

app.dimestcv_centeratdimpcaCheckBox.Value = est_params_cv.num_comp_center_around_dim_pca;
app.dimestcv_repsEditField.Value = est_params_cv.reps;
app.dimestcv_includeshuffversionCheckBox.Value = est_params_cv.include_shuff_version;

%% default ensemble params

app.ens_ensemethodDropDown.Value = ens_params.ensamble_method;
app.ens_normalizeDropDown.Value = ens_params.normalize;
app.ens_ensextractionDropDown.Value = ens_params.ensamble_extraction;
app.ens_extractionthreshDropDown.Value = ens_params.ensamble_extraction_thresh;
app.ens_signalzthreshEditField.Value = ens_params.signal_z_thresh;
app.ens_shuffthreshprcEditField.Value = ens_params.shuff_thresh_percent;
app.ens_hclustmethodDropDown.Value = ens_params.hcluster_method;
app.ens_hclustmetricDropDown.Value = ens_params.hcluster_distance_metric;
app.ens_plotstuffCheckBox.Value = ens_params.plot_stuff;
 
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