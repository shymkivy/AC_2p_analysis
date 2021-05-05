function f_dv_update_dset_info(app)

app.mplSpinner.Value = min([app.mplSpinner.Value, app.ddata.num_planes]);
app.MMNfreqEditField.Value = num2str(app.ddata.MMN_freq{1});

n_pl = app.mplSpinner.Value;

deconv_methods = {};
if ~isempty(app.ddata.OA_data{n_pl}.est.S)
    deconv_methods = [deconv_methods, {'OA_deconv'}];
end
if ~isempty(app.ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.S)
    deconv_methods = [deconv_methods, {'smooth_dfdt'}];
end
if ~isempty(cat(1,app.ddata.OA_data{n_pl}.proc.deconv.MCMC.S{:}))
    deconv_methods = [deconv_methods, {'MCMC'}];
end
if ~isempty(cat(1,app.ddata.OA_data{n_pl}.proc.deconv.c_foopsi.S{:}))
    deconv_methods = [deconv_methods, {'constrained_foopsi'}];
end

app.DeconvolutionmethodDropDown.Items = deconv_methods;

f_dv_update_A(app);
f_dv_update_cell(app);


end