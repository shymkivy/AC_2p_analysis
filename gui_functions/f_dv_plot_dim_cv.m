function f_dv_plot_dim_cv(app)

ddata = app.ddata;
data_dim_cv = ddata.data_dim_cv{1};

%%
if ~isempty(data_dim_cv)
    f_plot_cv_error_3D(data_dim_cv.est_params_list, data_dim_cv.est_params_list_shuff, 'smooth_SD', 'num_comp', 'test_err');
    ax1 = gca;
    ax1.Title.String = sprintf('dset %s', app.ddata.dset_name_full{1});          
else
    disp('Have not computed dim with CV yet...');
end

end