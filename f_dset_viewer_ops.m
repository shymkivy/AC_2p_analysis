function f_dset_viewer_ops(app)

app.gui_ops.reg_data_path = ...
    'C:\Users\ys2605\Desktop\stuff\AC_data\wf_registration_data\reg_save_6_10_21.mat';

app.gui_ops.mat_data_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\echo_save_12_19_21.mat.mat';


%% stat default params

stats.stat_method = 'shuff_pool'; % 'shuff_pool', 'shuff_locwise', 'z_thresh'
stats.stat_source = 'Freqs_dd'; % 'All', 'Freqs', 'Freqs_dd'
stats.z_thresh = 3;
stats.peak_bin_time = 0.250; % in sec
stats.num_shuff_samp = 2000;
stats.base_resp_win = [-2 3];
stats.loco_thresh = 99; % in percent;

%% ensemble default params




app.gui_ops.stats = stats;

end