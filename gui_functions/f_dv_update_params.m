function f_dv_update_params(app)

ddata = app.ddata;
n_pl = app.mplSpinner.Value;
stats1 = ddata.stats{n_pl};

if ~isempty(stats1)
    app.currentZthreshEditField.Value = stats1.stat_params.z_thresh;
    app.currentPvalEditField.Value = 1 - normcdf(app.currentZthreshEditField.Value);
    app.currentstatsmethodEditField.Value = stats1.stat_params.stat_method;
    app.currentstatssourceEditField.Value = stats1.stat_params.stat_source;
    app.currentPeakbintimesecEditField.Value = stats1.stat_params.peak_bin_time;
    app.currentNumshuffsampEditField.Value = stats1.stat_params.num_shuff_samp;
    app.currentbasewinEditField.Value = stats1.stat_params.base_resp_win(1);
    app.currentrespwinEditField.Value = stats1.stat_params.base_resp_win(2);
    app.currentLocothreshEditField.Value = stats1.stat_params.loco_thresh;
end

end