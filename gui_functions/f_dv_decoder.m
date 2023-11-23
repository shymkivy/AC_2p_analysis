function f_dv_decoder(app)

tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';

ddata = app.ddata;
cdata = f_dv_get_cdata(app);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[~, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1).volume_period);

firing_rate = cat(1,cdata.S_sm);

stats1 = cat(1,ddata.stats{:});

resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all);
resp_cells2 = logical(sum(resp_cells,2));

resp_cells3 = resp_cells(resp_cells2,:);

firing_rate2 = firing_rate(resp_cells2,:);

stim_times = ddata.stim_frame_index{1};
mmn_freq = ddata.MMN_freq{1};

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[~, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1).volume_period);

trial_types = app.ddata.trial_types{1};
trial_data_sort = f_get_stim_trig_resp(firing_rate2, stim_times, trial_frames);

trial_types2 = trial_types(2:400);

trial_data_sort2 = trial_data_sort(:,:,2:400);
trial_data_sort3 = squeeze(mean(trial_data_sort2,2));

trial_data_sort3n = trial_data_sort3 - mean(trial_data_sort3(:));
trial_data_sort3n = trial_data_sort3n/std(trial_data_sort3n(:));


[trainedClassifier_tree, validationAccuracy_tree] = decoder_tree(trial_data_sort3n', trial_types2);

[trainedClassifier_bayes, validationAccuracy_bayes] = decoder_naivebayes(trial_data_sort3n', trial_types2);

[trainedClassifier_svm, validationAccuracy_svm] = decoder_svm(trial_data_sort3n', trial_types2, 1);

trial_types2_shuff = trial_types2(randperm(numel(trial_types2)));

[trainedClassifier_tree_shuff, validationAccuracy_tree_shuff] = decoder_tree(trial_data_sort3n', trial_types2_shuff);

[trainedClassifier_bayes_shuff, validationAccuracy_bayes_shuff] = decoder_naivebayes(trial_data_sort3n', trial_types2_shuff);

[trainedClassifier_svm_shuff, validationAccuracy_svm_shuff] = decoder_svm(trial_data_sort3n', trial_types2_shuff, 1);

trial_types2_last = trial_types(1:399);

[trainedClassifier_tree_last, validationAccuracy_tree_last] = decoder_tree(trial_data_sort3n', trial_types2_last);

[trainedClassifier_bayes_last, validationAccuracy_bayes_last] = decoder_naivebayes(trial_data_sort3n', trial_types2_last);

[trainedClassifier_svm_last, validationAccuracy_svm_last] = decoder_svm(trial_data_sort3n', trial_types2_last, 1);



% mv regression
MnrModel_cont = fitmnr(trial_data_sort3n', trial_types2);


[d,p,stats] = manova1(trial_data_sort3n', trial_types2);

load fisheriris

MnrModel = fitmnr(meas,species);
MnrModel.Coefficients


% beta = mvregress(X,Y)

end