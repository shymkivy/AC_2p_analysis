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
params.z_thresh_new = app.ZthreshnewEditField.Value;
params.stat_source = app.StatsourceDropDown.Value;

end