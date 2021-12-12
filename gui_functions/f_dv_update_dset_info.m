function f_dv_update_dset_info(app)

idx1 = app.current_data_idx;

ddata = app.data(idx1,:);
%%
app.mplSpinner.Value = min([app.mplSpinner.Value, ddata.num_planes]);
app.MMNfreqEditField.Value = num2str(ddata.MMN_freq{1});
app.FramerateEditField.Value = 1000/ddata.proc_data{1}.frame_data.volume_period;

%%
frame_period = 1/app.FramerateEditField.Value;
trial_window = app.working_ops.trial_window;
trial_window_t = (ceil(trial_window(1)/frame_period):floor(trial_window(2)/frame_period))*frame_period;

app.working_ops.trial_window_t = trial_window_t;
app.working_ops.trial_num_baseline_resp_frames = [sum(trial_window_t<=0) sum(trial_window_t>0)];     

%%
n_pl = app.mplSpinner.Value;

deconv_methods = {};
if ~isempty(ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.S)
    deconv_methods = [deconv_methods, {'smooth_dfdt'}];
end
if ~isempty(cat(1,ddata.OA_data{n_pl}.proc.deconv.MCMC.S{:}))
    deconv_methods = [deconv_methods, {'MCMC'}];
end
if ~isempty(ddata.OA_data{n_pl}.est.S)
    deconv_methods = [deconv_methods, {'OA_deconv'}];
end
if ~isempty(cat(1,ddata.OA_data{n_pl}.proc.deconv.c_foopsi.S{:}))
    deconv_methods = [deconv_methods, {'constrained_foopsi'}];
end

app.DeconvolutionmethodDropDown.Items = deconv_methods;

%% gather C and S
params = f_dv_gather_params(app);
params.ddata = ddata;

% data_gr = app.SelectdatagroupButtonGroup.SelectedObject.Text;
% if strcmpi(data_gr, 'plane')
%     app.cdata = {f_dv_compute_cdata(app, params)};
% else
    cdata_all = cell(5,1);
    for n_pl2 = 1:ddata.num_planes
        params.n_pl = n_pl2;
        cdata_all{n_pl2} = f_dv_compute_cdata(app, params);
    end
    app.cdata = cat(1, cdata_all{:});
%end

%% compute dset statistics
for n_pl2 = 1:ddata.num_planes
    if isempty(ddata.stats{n_pl2})
        params.cdata = app.cdata(n_pl2,:);
        params.n_pl = n_pl2;
        app.data(idx1,:).stats{n_pl2} = f_dv_compute_stats_core(app, params);
        app.ddata = app.data(idx1,:);
    end
end

if isfield(app.data(idx1,:).stats{n_pl}, 'z_thresh')
    app.ZthreshcurrentEditField.Value = app.data(idx1,:).stats{n_pl}.z_thresh;
else 
    app.ZthreshcurrentEditField.Value = app.data(idx1,:).stats{n_pl}.stat_params.z_thresh;
end
app.pvalcurrentEditField.Value = 1 - normcdf(app.ZthreshcurrentEditField.Value);

if ~isempty(ddata.data_dim_pca{n_pl})
    app.DimpcaEditField.Value = ddata.data_dim_pca{n_pl}.dimensionality_corr;
else
    app.DimpcaEditField.Value = 0;
end

if ~isempty(ddata.data_dim_cv{n_pl})
    app.DimCVEditField.Value = ddata.data_dim_cv{n_pl}.dimensionality_corr;
else
    app.DimCVEditField.Value = 0;
end
end