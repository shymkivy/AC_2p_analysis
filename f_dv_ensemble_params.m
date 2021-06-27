function params = f_dv_ensemble_params(~, corr_dim)
%%
est_params_pca.normalize = 'norm_mean_std';
est_params_pca.dim_est_num_reps = 50;
est_params_pca.plot_stuff = 0;

%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles
% **params** are best params

est_params_cv.ensamble_method = 'pca';              % options: svd, pca (faster than svd), nmf, ica                % SVD is most optimal for encoding, NMF rotates components into something that is real and interpretable
est_params_cv.normalize = 'norm_mean_std'; % **'norm_mean_std'**, 'norm_mean' 'none'   % either way, need to normalize the power of signal in each cell, otherwise dimred will pull out individual cells
est_params_cv.shuffle_data_chunks = 1;   % 1 or 0, keeping cell correlations   % if the sequence of trial presentation contains information, you will need to shuffle. Also need to do in chunks because adjacent time bins are slightly correlated
% ---- input one or range of values to estimate across following
est_params_cv.smooth_SD = 0;       % larger window will capture 'sequences' of ensembles, if window is smaller than optimal, you will end up splitting those into more components
if exist('cv_corr_dim', 'var')
    est_params_cv.num_comp = (corr_dim-10):2:(corr_dim+10); 
else
    est_params_cv.num_comp = 1:2:30;
end
est_params_cv.reps = 5;              % how many repeats per param 
est_params_cv.include_shuff_version = 0;

%%
% NMF ensemble detection is best
% for NMF best to use norm_rms(keep values positive), otherwise can also use norm_mean_std
% NMF 14 comp
% SVD 11-14 comp?
ens_params.ensamble_method = 'nmf'; % options: svd, **nmf**, ica     % here NMF is
ens_params.num_comp = corr_dim;
ens_params.smooth_SD = 0; % 110 is better?
ens_params.normalize = 'norm_mean_std'; % 'norm_mean_std', 'norm_mean' 'none'
ens_params.ensamble_extraction = 'thresh'; %  **'thresh'(only for nmf)** 'clust'(for all)
% --- for thresh detection (only nmf)
ens_params.ensamble_extraction_thresh = 'signal_z'; % 'shuff' 'signal_z' 'signal_clust_thresh'
ens_params.signal_z_thresh = 2;
ens_params.shuff_thresh_percent = 95;
% --- for clust detection and general sorting 
ens_params.hcluster_method = 'average';  % ward(inner square), **average**, single(shortest)     
ens_params.hcluster_distance_metric = 'cosine';  % none, euclidean, squaredeuclidean, **cosine**, hammilarity, rbf% for low component number better euclidean, otherwise use cosine
ens_params.corr_cell_thresh_percent = 95;   % to remove cells with no significant correlations
% --- other
ens_params.plot_stuff = 0;
ens_params.acc_shuff_reps = 20;

%%
params.est_params_pca = est_params_pca;
params.est_params_cv = est_params_cv;
params.ens_params = ens_params;
end