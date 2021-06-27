function data_dim_cv = f_dv_estimate_dim_cv_core(params)

cv_corr_dim = round(params.data_dim_pca.dimensionality_corr);

%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles
params = f_dv_ensemble_params(app, cv_corr_dim);
est_params_cv = params.est_params_cv;
%%
firing_rate = f_dv_get_firing_rate(params.cdata);

%%
est_params_cv.n_rep = 1:est_params_cv.reps;
est_params_list = f_build_param_list(est_params_cv, {'smooth_SD', 'num_comp', 'n_rep'});
if est_params_cv.include_shuff_version
    est_params_list_s = est_params_list;
end

volume_period = params.cdata.volume_period;

%%
%firing_rate = f_dv_get_current_firing_rate(app);
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

    fig1 = f_plot_cv_error_3D(est_params_list, est_params_list_s, 'smooth_SD', 'num_comp', 'test_err');
else
    fig1 = f_plot_cv_error_3D(est_params_list, [], 'smooth_SD', 'num_comp', 'test_err');
end

fig1.Children(2).Title.String = sprintf('dset %s', params.ddata.experiment{1});          

data_dim_cv.test_err_data = test_err_data;
data_dim_cv.num_comp = est_params_cv.num_comp';
data_dim_cv.dimensionality_corr = est_params_cv.num_comp(min_ind);

end