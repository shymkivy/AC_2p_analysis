function params = f_dv_gather_params(app)

params.n_pl = app.mplSpinner.Value;

if strcmpi(app.SelectdatagroupButtonGroup.SelectedObject.Text, 'plane')
    params.planes = app.mplSpinner.Value;
else
    params.planes = 1:app.data(app.current_data_idx,:).num_planes;
end

%%
params.n_dset = find(app.current_data_idx);
params.deconvolution = app.DeconvolutionmethodDropDown.Value;
params.smooth = app.SmoothCheckBox.Value;
params.smooth_sigma = app.SmoothsigmamsEditField.Value;
params.rectify_spikes = app.RectifyspikesCheckBox.Value;
params.subtract_mean_spikes = app.SubtractmeanspikesCheckBox.Value;
params.normalize_max_spikes = app.NormalizemaxspikesCheckBox.Value;

%%
stats.stat_source = app.stats_StatsourceDropDown.Value;
stats.stat_method = app.stats_StatmethodDropDown.Value;
stats.z_thresh = app.stats_ZthreshEditField.Value;
stats.peak_bin_time = app.stats_PeakbintimesecEditField.Value;
stats.num_shuff_samp = app.stats_NumshuffsampEditField.Value;
stats.base_resp_win = f_str_to_array(app.stats_BaserespwindowEditField.Value);
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
ens_params.ensamble_method = app.ens_ensemethodDropDown.Value;
ens_params.normalize = app.ens_normalizeDropDown.Value;
ens_params.ensamble_extraction = app.ens_ensextractionDropDown.Value;
ens_params.ensamble_extraction_thresh = app.ens_extractionthreshDropDown.Value;
ens_params.signal_z_thresh = app.ens_signalzthreshEditField.Value;
ens_params.shuff_thresh_percent = app.ens_shuffthreshprcEditField.Value;
ens_params.hcluster_method = app.ens_hclustmethodDropDown.Value;
ens_params.hcluster_distance_metric = app.ens_hclustmetricDropDown.Value;
ens_params.plot_stuff = app.ens_plotstuffCheckBox.Value;

%%
params.stats = stats;
params.est_params_pca = est_params_pca;
params.est_params_cv = est_params_cv;
params.ens_params = ens_params;

end