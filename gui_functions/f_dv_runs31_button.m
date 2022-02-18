function f_dv_runs31_button(app)
%% load preprocessing parameters
ops = app.gui_ops.ops;

%% from S31 file s31MPL_audio_mmn_analysis.m
ops.processed_data_tag = 'processed_data';
ops.OA_output_tag = 'results_sort';

%% List of files to load

%which files to analyze
%AC_files_preOA_7_31_19
%AC_files_OA_10_27_19
%AC_files_MPL_files_12_11_19
%AC_files_OA_3_16_20;
%AC_files_OA_5_23_20_lab;

AC_files_xls_import;

%%
ops = f_mpl_process_ops(ops);

%%
[data, ops] = f_mpl_load_data(ops);

%%
data = f_dv_preprocess_data(data, ops);

app.data = data;
app.ops = ops;

%%
f_dv_compute_cdata_all(app)

%% additional 
f_dv_initialize_post_load(app);
end