function f_dv_ensemble_analysis(app)

ddata = app.ddata;

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
ens_params.num_comp = app.numenstofindEditField.Value;
ens_params.smooth_SD = app.ens_extrasmoothSDEditField.Value;

%%
ens_params.vol_period = ddata.proc_data{1}.frame_data.volume_period;

%%
cdata = f_dv_get_cdata(app);
firing_rate = cat(1,cdata.S_sm);
active_cells = sum(firing_rate,2) ~= 0;
firing_rate(~active_cells,:) = [];

% num_cells = size(firing_rate,1);
% firing_rate = firing_rate(randperm(num_cells),:);

%% extract ensembles
firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/ens_params.vol_period);
firing_rate_sm_norm = f_normalize(firing_rate_sm, ens_params.normalize);
ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm_norm, ens_params);

%% evaluate components
firing_rate_norm = f_normalize(firing_rate, ens_params.normalize);
acc_out_d = f_evaluate_ens_cv(ens_out, firing_rate_norm, ens_params);

ens_params.acc_out_d = acc_out_d;
ens_params.ens_out = ens_out;

ddata_idx = strcmpi(app.ddata.dset_name_full, app.data.dset_name_full);
app.data(ddata_idx,:).ensembles{1} = ens_params;
app.ddata.ensembles{1} = ens_params;

disp('Done');

end