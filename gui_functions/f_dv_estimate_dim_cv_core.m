function data_dim_cv = f_dv_estimate_dim_cv_core(params)

corr_dim = round(params.data_dim_pca.dimensionality_corr);

%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles
% **params** are best params
include_shuff_version = 0;
est_params.ensamble_method = 'pca';              % options: svd, pca (faster than svd), nmf, ica                % SVD is most optimal for encoding, NMF rotates components into something that is real and interpretable
est_params.normalize = 'norm_mean_std'; % **'norm_mean_std'**, 'norm_mean' 'none'   % either way, need to normalize the power of signal in each cell, otherwise dimred will pull out individual cells
est_params.shuffle_data_chunks = 1;   % 1 or 0, keeping cell correlations   % if the sequence of trial presentation contains information, you will need to shuffle. Also need to do in chunks because adjacent time bins are slightly correlated
% ---- input one or range of values to estimate across following
est_params.smooth_SD = 0;       % larger window will capture 'sequences' of ensembles, if window is smaller than optimal, you will end up splitting those into more components
est_params.num_comp = (corr_dim-5):1:(corr_dim+5);       
est_params.reps = 4;              % how many repeats per param 

%%
firing_rate = f_dv_get_firing_rate(params.cdata);

%%
est_params.n_rep = 1:est_params.reps;
est_params_list = f_build_param_list(est_params, {'smooth_SD', 'num_comp', 'n_rep'});
if include_shuff_version
    est_params_list_s = est_params_list;
end

volume_period = params.cdata.volume_period;

%%
%firing_rate = f_dv_get_current_firing_rate(app);
firing_rate_norm = f_normalize(firing_rate, est_params.normalize);

%% estimate num of components

est_params_list = f_ens_estimate_dim_params(firing_rate_norm, est_params_list, volume_period);
test_err_data = reshape([est_params_list.test_err]', [], est_params_list(1).reps);
[~, min_ind] = min(mean(test_err_data,2));
sd_all = mean(test_err_data,2);
fprintf('Optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', sd_all(min_ind), est_params.ensamble_method, est_params_list(min_ind).num_comp);

if include_shuff_version
    fprintf('Now shuff version...\n');
    firing_rate_norm_s = f_shuffle_data(firing_rate_norm);
    est_params_list_s = f_ens_estimate_dim_params(firing_rate_norm_s, est_params_list_s, volume_period);
    [~, min_ind] = min(mean(reshape([est_params_list_s.test_err]', [], est_params_list_s(1).reps),2));
    sd_all = mean(reshape([est_params_list_s.smooth_SD]', [], est_params_list_s(1).reps),2);
    fprintf('For shuff, optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', sd_all(min_ind), est_params.ensamble_method, est_params_list(min_ind).num_comp);

    fig1 = f_plot_cv_error_3D(est_params_list, est_params_list_s, 'smooth_SD', 'num_comp', 'test_err');
else
    fig1 = f_plot_cv_error_3D(est_params_list, [], 'smooth_SD', 'num_comp', 'test_err');
end

fig1.Children(2).Title.String = sprintf('dset %s', params.ddata.experiment{1});          

data_dim_cv.test_err_data = test_err_data;
data_dim_cv.num_comp = est_params.num_comp';
data_dim_cv.dimensionality_corr = est_params.num_comp(min_ind);

end