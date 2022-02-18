function f_dv_update_dset_info(app)

idx1 = app.current_data_idx;

ddata = app.data(idx1,:);
%%
app.mplSpinner.Value = min([app.mplSpinner.Value, ddata.num_planes]);
app.MMNfreqEditField.Value = num2str(ddata.MMN_freq{1});
app.FramerateEditField.Value = 1000/ddata.proc_data{1}.frame_data.volume_period;

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

cdata_all = cell(5,1);
for n_pl2 = 1:ddata.num_planes
    params.n_pl = n_pl2;
    cdata_all{n_pl2} = f_dv_compute_cdata(app, params);
end
app.cdata = cdata_all;

%% compute dset statistics
if isempty(ddata.stats{1})
    fprintf('Computing stats _/%d planes: ', ddata.num_planes);
end
for n_pl2 = 1:ddata.num_planes
    
    if isempty(ddata.stats{n_pl2})
        fprintf('%d..',n_pl2);
        params.cdata = cdata_all{n_pl2};
        params.n_pl = n_pl2;
        app.data(idx1,:).stats{n_pl2} = f_dv_compute_stats_core(app, params);
        app.ddata = app.data(idx1,:);
    end
end
if isempty(ddata.stats{1})
    fprintf(' Done\n');
end
if ~isempty(ddata.data_dim_pca{1})
    app.DimpcaEditField.Value = ddata.data_dim_pca{1}.dimensionality_corr;
else
    app.DimpcaEditField.Value = 0;
end

if ~isempty(ddata.data_dim_cv{1})
    app.DimCVEditField.Value = ddata.data_dim_cv{1}.dimensionality_corr;
else
    app.DimCVEditField.Value = 0;
end
end