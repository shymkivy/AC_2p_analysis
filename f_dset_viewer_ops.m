function f_dset_viewer_ops(app)

app.gui_ops.reg_data_path = ...
    'C:\Users\ys2605\Desktop\stuff\AC_data\wf_registration_data\reg_save_6_10_21.mat';

app.gui_ops.mat_data_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\save_6_27_21.mat';

app.gui_ops.save_var_list = {'stats', 'data_dim_pca', 'data_dim_cv', 'ensembles', 'ensemble_stats', 'ensemble_tuning'};

end