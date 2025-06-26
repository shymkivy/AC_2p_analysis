function params = f_dv_gather_params(app)

params.n_pl = app.mplSpinner.Value;

if ~isempty(app.data)
    current_dset = app.DatasetDropDown.Value;
    idx1 = strcmpi(app.data.dset_name_full, current_dset);
    params.current_data_idx = idx1;
    if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
        params.planes = app.mplSpinner.Value;
    else
        params.planes = 1:app.data(app.current_data_idx,:).num_planes;
    end
    params.n_dset = find(app.current_data_idx);
end

params.trial_type = app.trialtypeDropDown.Value;
params.region = app.regiontoplotDropDown.Value;
params.responsive_cells_type = app.ResposivecellstypeDropDown.Value;
params.responsive_cells_select = app.ResponsivecellsselectDropDown.Value;
params.responsive_thresh = app.RespthreshEditField.Value;
params.use_reg_data_labels = app.UseregdatalabelsCheckBox.Value;
params.pool_regions = app.poolregionsCheckBox.Value;
params.plot_feature = app.plotfeatureDropDown.Value;


params.convert_to_z = app.ConverttoZCheckBox.Value;
params.stats_between = app.statsbetweenDropDown.Value;
params.plot_stats = app.plotstatsCheckBox.Value;
params.plot_super_deets = app.plotsuperdeetsCheckBox.Value;
params.stim_window = app.winanalyzeDropDown.Value;
params.bar_plot_type = app.BarplottypeDropDown.Value;
params.marginalize_dist = app.MarginalizedistCheckBox.Value;
params.kde_smooth_factor = app.kdesmfactorEditField.Value;
params.plot_type = app.plottypeDropDown.Value;
params.max_y_lim = app.maxYlimEditField.Value;
params.min_y_lim = app.minYlimEditField.Value;
params.x_var = app.xvarDropDown.Value;
params.y_var = app.yvarDropDown.Value;
params.individual_tials = app.IndividualtrialsCheckBox.Value;

% plotting
params.plot_stim = app.PlotstimCheckBox.Value;
params.stim_transparancy = app.StimtranspEditField.Value;
params.stim_freq_color = app.stimcolorSpinner.Value;

if ~isempty(app.data)
    params.paradigm = app.data.paradigm{1};
    params.data_selection = app.SelectdatagroupDropDown.Value;
    params.current_dset_idx = strcmpi(app.data.dset_name_full, app.ddata.dset_name_full);
end

params.trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
params.onset_window = f_str_to_array(app.stats_OnsetRespwinEditField.Value);
params.offset_window = f_str_to_array(app.stats_OffsetRespwinEditField.Value);
params.plot_lims = f_str_to_array(app.plot_BaserespwinEditField.Value);

params.distance_method = lower(app.DistmethodDropDown.Value);
params.do_similarity = app.dosimilarityCheckBox.Value;
params.distance_reference = app.DistreferenceDropDown.Value;
params.normalize_euclidean = app.normalizeeucCheckBox.Value;
params.dim_red_method = app.DimredmethodDropDown.Value;
params.subtract_mean = app.subtractmeanCheckBox.Value;
params.scale_by_var = app.scalebyvarCheckBox.Value;
params.plot_subtracted_mean = app.plotsubmeanCheckBox.Value;
params.mat_tri = app.MattridataplotDropDown.Value;
params.colormap = app.ColormapDropDown.Value;
params.invert_cmap = app.InvertcmapCheckBox.Value;
params.non_handle_method = app.nanhandlemetDropDown.Value;

params.select_resp_cells = app.selectrespcellsCheckBox.Value;
params.resort_by_ens = app.resortbyensCheckBox.Value;
params.sort_trials = app.sorttrialsCheckBox.Value;
params.sort_with_full_firing_rate = app.sortwithfullfrCheckBox.Value;
params.shuffle_cells = app.shufflecellsCheckBox.Value;
params.shuffle_trials = app.shuffletrialsCheckBox.Value;
params.sort_by_trial_type = app.sortbytrialtypeCheckBox.Value;
params.sort_by_prev_trial = app.sortprevtrialCheckBox.Value;
params.sort_by_similarity = app.sortbycellsimilarityCheckBox.Value;
params.plot_sorting_stuff = app.plotsortingstuffCheckBox.Value;
params.sort_cell_by_trace_ave = app.sortcellbytraceaveCheckBox.Value;
params.plot_prev_trial = app.plotprevtrialCheckBox.Value;
params.plot_pop_vector = app.plotpopvecCheckBox.Value;
params.sort_ens_trial = app.SortenstrialsCheckBox.Value;

params.pca_distances = app.PCAdistancesCheckBox.Value;
params.plot_pca_dim = app.plotPCAdimSpinner.Value;

params.num_axes_plot = app.NumplotaxesDropDown.Value;
params.num_comp_plot = app.numcompplotSpinner.Value;

params.ld_control_pad = app.LDcontpadEditField.Value;

params.render_painters = app.FigrenderpaintersCheckBox.Value;
params.shadow_on3d = app.shadowon3dCheckBox.Value;
params.shadow_axis_locs = [app.FlipshadowXCheckBox.Value, app.FlipshadowYCheckBox.Value, app.FlipshadowZCheckBox.Value] + 1;
params.reverse_xyz = [app.ReverseXCheckBox.Value, app.ReverseYCheckBox.Value, app.ReverseZCheckBox.Value];
params.grid_on = app.gridon3dCheckBox.Value;



% ensless corr
params.samp_range_min = app.samprangeminEditField.Value;
params.samp_range_max = app.samprangemaxEditField.Value;
params.samp_range_interval = app.samprangeintervEditField.Value;
params.num_samples = app.sampnumEditField.Value;
params.equalize_samp_size = app.EqualizetrialsampsizeCheckBox.Value;

params.sort_type = app.SorttypeDropDown.Value; %'trials', 'cells'
params.corr_type = app.CorrtypeDropDown.Value; % 'SI_cosine', 'SI_correlation', 'pca_dim'

params.plot_ens = app.plotensamblesCheckBox.Value;
params.plot_fractions = app.plotfractionsCheckBox.Value;

% widefield
params.freq_im_num = app.numfreqSpinner.Value;
params.anchor_dset = app.anchordsetSpinner.Value;
params.contours = app.ContoursDropDown.Value;
params.size_factor = app.SizefactorEditField.Value;
params.plot_borders = app.AreabordersCheckBox.Value;
params.plot_nontuned = app.PlotnontunedCheckBox.Value;
params.white_bkg = app.WhitebkgCheckBox.Value;
params.new_plots = app.NewplotsCheckBox.Value;
params.size_mag_adjust = app.SizemagadjustCheckBox.Value;
%%
params.deconvolution = app.DeconvolutionmethodDropDown.Value;
params.smooth = app.SmoothCheckBox.Value;
params.smooth_sigma = app.SmoothsigmamsEditField.Value;
params.rectify_spikes = app.RectifyspikesCheckBox.Value;
params.subtract_mean_spikes = app.SubtractmeanspikesCheckBox.Value;
params.normalize_max_spikes = app.NormalizemaxspikesCheckBox.Value;

params.anchor_reg_dset = app.anchordsetSpinner.Value;


%%
stats.stat_source = app.stats_StatsourceDropDown.Value;
stats.stat_method = app.stats_StatmethodDropDown.Value;
stats.z_thresh = app.stats_ZthreshEditField.Value;
stats.peak_bin_time = app.stats_PeakbintimesecEditField.Value;
stats.num_shuff_samp = app.stats_NumshuffsampEditField.Value;
stats.base_resp_win = f_str_to_array(app.stats_BaserespwinEditField.Value);
stats.lim_sig_resp_win = f_str_to_array(app.stats_LimSigRespwinEditField.Value);
stats.resp_win_onset = f_str_to_array(app.stats_OnsetRespwinEditField.Value);
stats.resp_win_offset = f_str_to_array(app.stats_OffsetRespwinEditField.Value);
stats.resp_win_middle = f_str_to_array(app.stats_MiddleRespwinEditField.Value);
stats.loco_thresh = app.stats_LocothreshEditField.Value;

%%
est_params_pca.normalize = app.dimestpca_normalizeDropDown.Value;
est_params_pca.shuffle_method = app.dimestpca_shufflemethodDropDown.Value;
est_params_pca.dim_est_num_reps = app.dimestpca_dimestnumrepsEditField.Value;
est_params_pca.total_dim_thresh = app.dimestpca_totaldimthreshEditField.Value;
est_params_pca.plot_stuff = app.dimestpca_PlotstuffCheckBox.Value;

%%
est_params_cv.normalize = app.dimestcv_normalizeDropDown.Value;
est_params_cv.ensemble_method = app.dimestcv_ensmethodDropDown.Value;
est_params_cv.shuffle_data_chunks = app.dimestcv_shuffledatachunksCheckBox.Value;
est_params_cv.smooth_SD = f_dv_get_dim_est_cv_sm(app);
est_params_cv.num_comp = f_dv_get_dim_est_cv_nc(app);
est_params_cv.center_at_dim_pca = app.dimestcv_centeratdimpcaCheckBox.Value;
est_params_cv.reps = app.dimestcv_repsEditField.Value;
est_params_cv.include_shuff_version = app.dimestcv_includeshuffversionCheckBox.Value;

%%
ens_params.ensemble_method = app.ens_ensemethodDropDown.Value;
ens_params.normalize = app.ens_normalizeDropDown.Value;
ens_params.ensemble_extraction = app.ens_ensextractionDropDown.Value;
ens_params.ensemble_extraction_thresh = app.ens_extractionthreshDropDown.Value;
ens_params.signal_z_thresh = app.ens_signalzthreshEditField.Value;
ens_params.shuff_thresh_percent = app.ens_shuffthreshprcEditField.Value;
ens_params.hcluster_method = app.ens_hclustmethodDropDown.Value;
ens_params.hcluster_distance_metric = app.ens_hclustmetricDropDown.Value;
ens_params.plot_stuff = app.ens_plotstuffCheckBox.Value;
ens_params.acc_shuff_reps = app.ensshuffrepsEditField.Value;
ens_params.smooth_SD = app.ens_extrasmoothSDEditField.Value;

%%
params.stats = stats;
params.est_params_pca = est_params_pca;
params.est_params_cv = est_params_cv;
params.ens_params = ens_params;

end