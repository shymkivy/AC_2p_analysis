function f_dset_viewer_ops(app)

app.gui_ops.reg_data_path = ...
    'C:\Users\ys2605\Desktop\stuff\AC_data\wf_registration_data\reg_save_6_10_21.mat';

app.gui_ops.mat_data_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\echo_save_12_19_21.mat.mat';


%%
app.gui_ops.stat_source_default = 'Freqs_dd'; % 'All', 'Freqs', 'Freqs_dd'

end