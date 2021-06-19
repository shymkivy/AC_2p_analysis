%% analyzing ensembles with NMF
% load data into 'firing_rate' variable (cells x frames)
clear;
close all;
addpath([pwd '\functions'])

%firing_rate = put data here

frame_rate = 30;

%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles
% **params** are best params
estimate_params = 1;    % do estimation?
include_shuff_version = 1;
est_params.ensamble_method = 'svd';              % options: svd, nmf, ica                % SVD is most optimal for encoding, NMF rotates components into something that is real and interpretable
est_params.normalize = 'norm_mean_std'; % **'norm_mean_std'**, 'norm_mean' 'none'   % either way, need to normalize the power of signal in each cell, otherwise dimred will pull out individual cells
est_params.shuffle_data_chunks = 0;   % 1 or 0, keeping cell correlations   % if the sequence of trial presentation contains information, you will need to shuffle. Also need to do in chunks because adjacent time bins are slightly correlated
% ---- input one or range of values to estimate across following
est_params.smooth_SD = 100;       % larger window will capture 'sequences' of ensembles, if window is smaller than optimal, you will end up splitting those into more components
est_params.num_comp = 1:2:20;       
est_params.reps = 2;              % how many repeats per param 

%%
est_params.n_rep = 1:est_params.reps;
est_params_list = f_build_param_list(est_params, {'smooth_SD', 'num_comp', 'n_rep'});
if include_shuff_version
    est_params_list_s = est_params_list;
end

%% input paramseters for ensemble analysis
% NMF ensemble detection is best
% for NMF best to use norm_rms(keep values positive), otherwise can also use norm_mean_std
% NMF 14 comp
% SVD 11-14 comp?
ens_params.ensamble_method = 'nmf'; % options: svd, **nmf**, ica     % here NMF is
ens_params.num_comp = 11;
ens_params.smooth_SD = 100; % 110 is better?
ens_params.normalize = 'norm_mean_std'; % 'norm_mean_std', 'norm_mean' 'none'
ens_params.ensamble_extraction = 'thresh'; %  **'thresh'(only for nmf)** 'clust'(for all)
ens_params.ensamble_extraction_thresh = 'shuff'; % 'shuff' 'signal_z' 'signal_clust_thresh'
% --- for thresh detection (only nmf)
ens_params.ensamble_extraction_thresh = 'signal_z'; % 'shuff' 'signal_z' 'signal_clust_thresh'
ens_params.signal_z_thresh = 2.5;
ens_params.shuff_thresh_percent = 95;
% --- for clust detection and general sorting 
ens_params.hcluster_method = 'average';  % ward(inner square), **average**, single(shortest)     
ens_params.hcluster_distance_metric = 'cosine';  % none, euclidean, squaredeuclidean, **cosine**, hammilarity, rbf% for low component number better euclidean, otherwise use cosine
ens_params.corr_cell_thresh_percent = 95;   % to remove cells with no significant correlations
% --- other
ens_params.plot_stuff = 1;

%%
vol_period = 1000/frame_rate;

%% remove inactive cells

active_cells = sum(firing_rate,2) > 0;
firing_rate(~active_cells,:) = [];

num_cells = size(firing_rate,1);

firing_rate = firing_rate(randperm(num_cells),:);
%% estimate best smoothing window
%% estimate smooth
if estimate_params
    fprintf('Estimating params n/%d reps: ',numel(est_params_list));
    %dim_corr = zeros(numel(estimate_smooth_list),1);
    for n_par = 1:numel(est_params_list)
        fprintf('%d..',n_par);

        firing_rate_norm = f_normalize(firing_rate, est_params.normalize);

        params1 = est_params_list(n_par);
        params1.vol_period = vol_period;
        accuracy = f_ens_estimate_corr_dim_cv(firing_rate_norm, params1);

        temp_fields = fields(accuracy);
        for n_fl = 1:numel(temp_fields)
            est_params_list(n_par).(temp_fields{n_fl}) = accuracy.(temp_fields{n_fl});
        end
    end
    fprintf('\nDone\n');
    [~, min_ind] = min([est_params_list.test_err]);
    fprintf('From provided range, optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', est_params_list(min_ind).smooth_SD, est_params.method, est_params_list(min_ind).num_comp);

    f_plot_cv_error_3D(est_params_list, [], 'smooth_SD', 'num_comp', 'test_err');
    ax1 = gca;
    ax1.Title.String = sprintf('%s, dset%d; %s L2 error from raw, (%s)',cond_name,n_dset,est_params.ensamble_method, ax1.Title.String);          
end

%% Smooth data
firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/vol_period);

%% extract ensambles
ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm, ens_params);