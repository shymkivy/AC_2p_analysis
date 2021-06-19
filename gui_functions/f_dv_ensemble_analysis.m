function f_dv_ensemble_analysis(app)

%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles
% **params** are best params
estimate_params = 1;    % do estimation?
include_shuff_version = 1;
est_params.ensamble_method = 'pca';              % options: svd, pca (faster than svd), nmf, ica                % SVD is most optimal for encoding, NMF rotates components into something that is real and interpretable
est_params.normalize = 'norm_mean_std'; % **'norm_mean_std'**, 'norm_mean' 'none'   % either way, need to normalize the power of signal in each cell, otherwise dimred will pull out individual cells
est_params.shuffle_data_chunks = 0;   % 1 or 0, keeping cell correlations   % if the sequence of trial presentation contains information, you will need to shuffle. Also need to do in chunks because adjacent time bins are slightly correlated
% ---- input one or range of values to estimate across following
est_params.smooth_SD = 100;       % larger window will capture 'sequences' of ensembles, if window is smaller than optimal, you will end up splitting those into more components
est_params.num_comp = 1:5:35;       
est_params.reps =1;              % how many repeats per param 

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
ens_params.num_comp = 21;
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
volume_period = app.ddata.proc_data{1}.frame_data.volume_period;
%%
n_pl = app.mplSpinner.Value;
tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';

stim_times = app.ddata.stim_frame_index{n_pl};
mmn_freq = app.ddata.MMN_freq{1};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
trial_types = app.ddata.trial_types{1};

firing_rate = app.cdata.S;

%%
active_cells = sum(firing_rate,2) ~= 0;
firing_rate(~active_cells,:) = [];

num_cells = size(firing_rate,1);

firing_rate = firing_rate(randperm(num_cells),:);

firing_rate_norm = f_normalize(firing_rate, est_params.normalize);
firing_rate_norm_s = f_shuffle_data(firing_rate_norm);

%% estimate best smoothing window
%% estimate smooth
if estimate_params

    est_params_list = f_ens_estimate_dim_params(firing_rate_norm, est_params_list, volume_period);
    [~, min_ind] = min(mean(reshape([est_params_list.test_err]', [], est_params_list(1).reps),2));
    sd_all = mean(reshape([est_params_list.smooth_SD]', [], est_params_list(1).reps),2);
    fprintf('Optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', sd_all(min_ind), est_params.ensamble_method, est_params_list(min_ind).num_comp);

    if include_shuff_version
        fprintf('Now shuff version...\n');
        est_params_list_s = f_ens_estimate_dim_params(firing_rate_norm_s, est_params_list_s, volume_period);
        [~, min_ind] = min(mean(reshape([est_params_list_s.test_err]', [], est_params_list_s(1).reps),2));
        sd_all = mean(reshape([est_params_list_s.smooth_SD]', [], est_params_list_s(1).reps),2);
        fprintf('For shuff, optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', sd_all(min_ind), est_params.ensamble_method, est_params_list(min_ind).num_comp);
        
        f_plot_cv_error_3D(est_params_list, est_params_list_s, 'smooth_SD', 'num_comp', 'test_err');
    else
        f_plot_cv_error_3D(est_params_list, [], 'smooth_SD', 'num_comp', 'test_err');
    end
    
    
    ax1 = gca;
    ax1.Title.String = sprintf('dset %s', app.ddata.experiment{1});          
end

%% Smooth data
firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/volume_period);
firing_rate_sm_s = f_shuffle_data(firing_rate_sm);

%% extract ensambles
ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm, ens_params);

%% evaluate components
ens_params.vol_period = volume_period;
acc_out_d = f_evaluate_ens_cv(ens_out, firing_rate_norm, ens_params);

num_shuff = 10;
ens_params_s = ens_params;
ens_params_s.num_comp = 2;
acc_out = cell(num_shuff,1);
for n_shuff = 1:num_shuff
    firing_rate_s = f_shuffle_data(firing_rate_norm);
    firing_rate_sm_s = f_smooth_gauss(firing_rate_s, ens_params_s.smooth_SD/vol_period);
    ens_out_s = f_ensemble_analysis_YS_raster(firing_rate_sm_s, ens_params_s);

    acc_out{n_shuff} = f_evaluate_ens_cv(ens_out_s, firing_rate_s, ens_params_s);
end



acc_out_full = cat(2,acc_out{:});
figure; plot(acc_out_full')

prctile(acc_out_full,1)
min(acc_out_full)

figure; hold on;
plot(acc_out_d)
shadedErrorBar_YS(1:ens_params_s.num_comp, nanmean(acc_out_full,2),nanstd(acc_out_full,[],2)/sqrt(num_shuff-1));

figure; plot(acc_out{1})

%% analyze ensembles
f_plot_raster_mean(firing_rate_sm(ens_out.ord_cell,:), 1);
title('raster cell sorted');

for n_comp = 1:numel(ens_out.cells.ens_list)
    cells1 = ens_out.cells.ens_list{n_comp};
    trials1 = ens_out.trials.ens_list{n_comp};
    scores1 = ens_out.cells.ens_scores(n_comp,:);

    f_plot_ensamble_deets(firing_rate_sm, cells1, trials1, scores1);
    title([ens_params.ensamble_method ' ensamble ' num2str(n_comp) '; ' num2str(tn')]);
end


%%

scores_3d = reshape(ens_out.cells.ens_scores, ens_params.num_comp, [], num_trials);

for n_comp = 1:ens_params.num_comp
    figure;
    for n_tr = 1:numel(tn)
        curr_tr = tn(n_tr);
        subplot(2,3,n_tr); hold on; axis tight
        traces1 = squeeze(scores_3d(n_comp,:,logical(trial_types_cut == curr_tr)));
        plot(traces1, 'color', [.5 .5 .5]);
        plot(mean(traces1,2), 'r', 'LineWidth', 2);
        ylim([0 0.22])
        if n_tr == 1
            title(['comp ' num2str(n_comp) '; ' num2str(tn')])
        end
    end
end


%% 
[~, max_ind] = max(firing_rate,[],2);
[~, seq1] = sort(max_ind, 'ascend');
figure; imagesc(firing_rate(seq1,:))

firing_rate_sh = f_shuffle_data(firing_rate);
[~, max_ind_s] = max(firing_rate_sh,[],2);
[~, seq1_s] = sort(max_ind_s, 'ascend');
figure; imagesc(firing_rate_sh(seq1_s,:))
title('shuffle')

end