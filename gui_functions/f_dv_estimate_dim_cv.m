function f_dv_estimate_dim_cv(app, params)

n_pl = params.n_pl;

if isempty(app.ddata.data_dim_est{n_pl})
    f_dv_estimate_dim_pca(app);
end

cv_corr_dim = round(app.ddata.data_dim_est{n_pl}.dimensionality_corr);

%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles
params = f_dv_ensemble_params(app, cv_corr_dim);
est_params_cv = params.est_params_cv;
%%
est_params_cv.n_rep = 1:est_params_cv.reps;
est_params_list = f_build_param_list(est_params_cv, {'smooth_SD', 'num_comp', 'n_rep'});
if est_params_cv.include_shuff_version
    est_params_list_s = est_params_list;
end

volume_period = app.ddata.proc_data{1}.frame_data.volume_period;
%%
% tn_all = f_dv_get_trial_number(app);
% tt_all = app.ops.context_types_all(tn_all)';
% 
% stim_times = app.ddata.stim_frame_index{n_pl};
% mmn_freq = app.ddata.MMN_freq{1};
% trig_window = app.working_ops.trial_num_baseline_resp_frames;
% trial_types = app.ddata.trial_types{1};


%%
firing_rate = f_dv_get_current_firing_rate(app);
firing_rate_norm = f_normalize(firing_rate, est_params_cv.normalize);

%% estimate num of components

est_params_list = f_ens_estimate_dim_params(firing_rate_norm, est_params_list, volume_period);
test_err_data = reshape([est_params_list.test_err]', [], est_params_list(1).reps);
[~, min_ind] = min(mean(test_err_data,2));
sd_all = mean(test_err_data,2);
fprintf('Optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', sd_all(min_ind), est_params_cv.ensamble_method, est_params_list(min_ind).num_comp);

if est_params_cv.include_shuff_version
    fprintf('Now shuff version...\n');
    firing_rate_norm_s = f_shuffle_data(firing_rate_norm);
    est_params_list_s = f_ens_estimate_dim_params(firing_rate_norm_s, est_params_list_s, volume_period);
    [~, min_ind] = min(mean(reshape([est_params_list_s.test_err]', [], est_params_list_s(1).reps),2));
    sd_all = mean(reshape([est_params_list_s.smooth_SD]', [], est_params_list_s(1).reps),2);
    fprintf('For shuff, optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', sd_all(min_ind), est_params_cv.ensamble_method, est_params_list(min_ind).num_comp);

    f_plot_cv_error_3D(est_params_list, est_params_list_s, 'smooth_SD', 'num_comp', 'test_err');
else
    f_plot_cv_error_3D(est_params_list, [], 'smooth_SD', 'num_comp', 'test_err');
end

ax1 = gca;
ax1.Title.String = sprintf('dset %s', app.ddata.experiment{1});          

data_dim_cv.test_err_data = test_err_data;
data_dim_cv.num_comp = est_params_cv.num_comp';
data_dim_cv.dimensionality_corr = est_params_cv.num_comp(min_ind);

%%
app.DimCVEditField.Value = data_dim_cv.dimensionality_corr;

max_planes = max(app.data.num_planes);
if ~sum(strcmpi(app.data.Properties.VariableNames, 'data_dim_cv'))
    app.data.data_dim_cv = cell(size(app.data,1),max_planes);
    app.ddata.data_dim_cv = cell(1,max_planes);
end

ddata_idx = strcmpi(app.ddata.experiment, app.data.experiment);
app.data(ddata_idx,:).data_dim_cv{n_pl} = data_dim_cv;
app.ddata.data_dim_cv{n_pl} = data_dim_cv;

disp('Done')

end