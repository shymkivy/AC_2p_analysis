function params = f_dv_gather_params(app)

params.n_pl = app.mplSpinner.Value;

if strcmpi(app.SelectdatagroupButtonGroup.SelectedObject.Text, 'plane')
    params.planes = app.mplSpinner.Value;
else
    params.planes = 1:app.data(app.current_data_idx,:).num_planes;
end

params.n_dset = find(app.current_data_idx);
params.deconvolution = app.DeconvolutionmethodDropDown.Value;
params.smooth = app.SmoothCheckBox.Value;
params.smooth_sigma = app.SmoothsigmamsEditField.Value;
params.rectify_spikes = app.RectifyspikesCheckBox.Value;
params.subtract_mean_spikes = app.SubtractmeanspikesCheckBox.Value;
params.normalize_max_spikes = app.NormalizemaxspikesCheckBox.Value;

params.stats.stat_source = app.newStatsourceDropDown.Value;
params.stats.stat_method = app.newStatmethodDropDown.Value;
params.stats.z_thresh = app.newZthreshEditField.Value;
params.stats.peak_bin_time = app.newPeakbintimesecEditField.Value;
params.stats.num_shuff_samp = app.newNumshuffsampEditField.Value;
params.stats.base_resp_win = [app.newbasewinEditField.Value app.newrespwinEditField.Value];
params.stats.loco_thresh = app.newLocothreshEditField.Value;

end