function ens_params = f_dv_ensemble_extract_core(~, params)

corr_dim = round(params.data_dim_pca.dimensionality_corr);

%% input paramseters for ensemble analysis
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
volume_period = params.cdata.volume_period;
%%
%n_pl = app.mplSpinner.Value;
%tn_all = f_dv_get_trial_number(app);
% tt_all = app.ops.context_types_all(tn_all)';
% 
% stim_times = app.ddata.stim_frame_index{n_pl};
% mmn_freq = app.ddata.MMN_freq{1};
% trig_window = app.working_ops.trial_num_baseline_resp_frames;
% trial_types = app.ddata.trial_types{1};

firing_rate = params.cdata.S;

%%
active_cells = sum(firing_rate,2) ~= 0;
firing_rate(~active_cells,:) = [];

%num_cells = size(firing_rate,1);
%firing_rate = firing_rate(randperm(num_cells),:);

firing_rate_norm = f_normalize(firing_rate, ens_params.normalize);
firing_rate_norm_s = f_shuffle_data(firing_rate_norm);

%% Smooth data
firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/volume_period);

%% extract ensambles
ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm, ens_params);

%% evaluate components
ens_params.vol_period = volume_period;
acc_out_d = f_evaluate_ens_cv(ens_out, firing_rate_norm, ens_params);

ens_params_s = ens_params;
ens_params_s.num_comp = 10;
ens_params_s.plot_stuff = 0;
acc_out_s = cell(ens_params.acc_shuff_reps,1);
fprintf('Computing shuff ens n/%d', ens_params.acc_shuff_reps);
for n_shuff = 1:ens_params.acc_shuff_reps
    firing_rate_sm_s = f_shuffle_data(firing_rate_sm);
    ens_out_s = f_ensemble_analysis_YS_raster(firing_rate_sm_s, ens_params_s);

    acc_out_s{n_shuff} = f_evaluate_ens_cv(ens_out_s, firing_rate_norm_s, ens_params_s);
    fprintf('--%d', n_shuff);
end
fprintf('\nDone\n');

try
    acc_out_full = cat(1,acc_out_s{:});
catch
    acc_out_full = [];
    for ii = 1:numel(acc_out_s)
        acc_out_full = [acc_out_full; acc_out_s{ii}(:)];
    end
    disp('error in acc extract');
end
    
thresh = prctile(acc_out_full(:), 1);

ens_params.acc_out_d = acc_out_d;
ens_params.acc_out_shuff = acc_out_s;
ens_params.acc_out_thresh = thresh;
ens_params.accepted_ensembles = acc_out_d<thresh;

% 
% figure; 
% subplot(2,1,1); hold on;
% [f_d,x_d] = ecdf(acc_out_d);
% [f_s,x_s] = ecdf(acc_out_full(:));
% plot(x_d, f_d, 'LineWidth', 2);
% plot(x_s, f_s, 'LineWidth', 2);
% line([thresh thresh], [0 1], 'Color', 'r', 'LineStyle', '--')
% legend('Data', 'Shuff', '99% Thresh')
% xlim([min(x_d) max(x_d)*2])
% bndw = 1;
% xlabel('test error');
% ylabel('component fraction')
% title('component test error ecdf')
% subplot(2,1,2); hold on;
% [f_d,x_d] = ksdensity(acc_out_d, 'Bandwidth', bndw);
% [f_s,x_s] = ksdensity(acc_out_full(acc_out_full<(max(x_d)*2)), 'Bandwidth', bndw);
% plot(x_d, f_d, 'LineWidth', 2);
% plot(x_s, f_s, 'LineWidth', 2);
% line([thresh thresh], [0 max([f_d, f_s])], 'Color', 'r', 'LineStyle', '--')
% legend('Data', 'Shuff', '99% Thresh')
% xlim([min(x_d) max(x_d)*2])
% ylim([0 max([f_d, f_s])])
% title('component test error kernel density')
% xlabel('test error');
% ylabel('component number')
% 
% fprintf('%.1f%% ensembles are above shuff\n', sum(x_d<thresh)/numel(x_d)*100);
% 
% %% analyze ensembles
% f_plot_raster_mean(firing_rate_sm(ens_out.ord_cell,:), 1);
% title('raster cell sorted');
% 
% for n_comp = 1:numel(ens_out.cells.ens_list)
%     cells1 = ens_out.cells.ens_list{n_comp};
%     trials1 = ens_out.trials.ens_list{n_comp};
%     scores1 = ens_out.cells.ens_scores(n_comp,:);
%     coeffs1 = ens_out.coeffs(cells1,n_comp);
% 
%     f_plot_ensamble_deets(firing_rate_sm, cells1, trials1, scores1, coeffs1);
%     title(sprintf('%s ensamble %d; test error=%.1f', ens_params.ensamble_method, n_comp, acc_out_d(n_comp)));
% end
% 
% 
% %%
% 
% figure; imagesc(ens_out.cells.ens_scores)
% ens_scores_sort = f_get_stim_trig_resp(ens_out.cells.ens_scores, stim_times, trig_window);
% [ens_scores_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(ens_scores_sort, trial_types, mmn_freq, app.ops);
% 
% % freqs
% for n_comp = 1:ens_params.num_comp
%     figure;
%     for n_tr = 1:10
%         curr_tt = app.ops.context_types_all(n_tr);
%         subplot(2,5,n_tr); hold on; axis tight
%         tr_idx = trial_types_wctx == curr_tt;
%         tr_ave = squeeze(ens_scores_sort_wctx(n_comp,:,tr_idx));
%         trace_mean = mean(tr_ave,2);
%         trace_sem = std(tr_ave,[],2)/sqrt(sum(tr_idx)-1);
%         plot(tr_ave, 'color', [.5 .5 .5]);
%         plot(trace_mean, 'r', 'LineWidth', 2);
%         ylim([0 0.22])
%         if n_tr == 1
%             title(['comp ' num2str(n_comp)])
%         end
%     end
% end
% 
% % ctx
% ctx_trials = [18 19 20 28 29 30];
% for n_comp = 1:ens_params.num_comp
%     figure;
%     for n_tr = 1:numel(ctx_trials)
%         curr_tt = app.ops.context_types_all(ctx_trials(n_tr));
%         subplot(2,3,n_tr); hold on; axis tight
%         tr_idx = trial_types_wctx == curr_tt;
%         tr_ave = squeeze(ens_scores_sort_wctx(n_comp,:,tr_idx));
%         trace_mean = mean(tr_ave,2);
%         trace_sem = std(tr_ave,[],2)/sqrt(sum(tr_idx)-1);
%         plot(tr_ave, 'color', [.5 .5 .5]);
%         plot(trace_mean, 'color', app.ops.context_types_all_colors2{ctx_trials(n_tr)}, 'LineWidth', 2);
%         ylim([0 0.22])
%         if n_tr == 1
%             title(['comp ' num2str(n_comp)])
%         end
%     end
% end


end