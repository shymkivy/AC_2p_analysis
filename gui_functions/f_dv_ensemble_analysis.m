function f_dv_ensemble_analysis(app)

ddata = app.ddata;
n_pl = app.mplSpinner.Value;

%%
corr_dim = round(ddata.data_dim_pca{n_pl}.dimensionality_corr);
params = f_dv_ensemble_params(app, corr_dim);

est_params = params.est_params_cv;
ens_params = params.ens_params;
%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles
estimate_params = 0;    % do estimation?

%%
est_params.n_rep = 1:est_params.reps;
est_params_list = f_build_param_list(est_params, {'smooth_SD', 'num_comp', 'n_rep'});
if est_params.include_shuff_version
    est_params_list_s = est_params_list;
end

%%
ens_params.vol_period = ddata.proc_data{1}.frame_data.volume_period;

%%

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

% num_cells = size(firing_rate,1);
% firing_rate = firing_rate(randperm(num_cells),:);

%% estimate best smoothing window
%% estimate smooth
if estimate_params
    firing_rate_norm = f_normalize(firing_rate, ens_params.normalize);
    est_params_list = f_ens_estimate_dim_params(firing_rate_norm, est_params_list, volume_period);
    [~, min_ind] = min(mean(reshape([est_params_list.test_err]', [], est_params_list(1).reps),2));
    sd_all = mean(reshape([est_params_list.smooth_SD]', [], est_params_list(1).reps),2);
    fprintf('Optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', sd_all(min_ind), est_params.ensamble_method, est_params_list(min_ind).num_comp);

    if est_params.include_shuff_version
        fprintf('Now shuff version...\n');
        est_params_list_s = f_ens_estimate_dim_params(firing_rate_sm_norm, est_params_list_s, volume_period);
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

%% extract ensambles
firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/ens_params.vol_period);
firing_rate_sm_norm = f_normalize(firing_rate_sm, ens_params.normalize);
ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm_norm, ens_params);

%% evaluate components
firing_rate_norm = f_normalize(firing_rate, ens_params.normalize);
acc_out_d = f_evaluate_ens_cv(ens_out, firing_rate_norm, ens_params);

%% shuffled data
ens_params.acc_shuff_reps = 5;
ens_params_s = ens_params;
ens_params_s.num_comp = 10;
ens_params_s.plot_stuff = 0;
acc_out_s = cell(ens_params.acc_shuff_reps,1);
fprintf('Computing shuff ens n/%d', ens_params.acc_shuff_reps);
for n_shuff = 1:ens_params.acc_shuff_reps
    firing_rate_s = f_shuffle_data(firing_rate);
    
    % estimation
    firing_rate_s_sm = f_smooth_gauss(firing_rate_s, ens_params.smooth_SD/ens_params.vol_period);
    firing_rate_s_sm_norm = f_normalize(firing_rate_s_sm, ens_params_s.normalize);
    ens_out_s = f_ensemble_analysis_YS_raster(firing_rate_s_sm_norm, ens_params_s);
    
    % accuracy
    firing_rate_s_norm = f_normalize(firing_rate_s, ens_params_s.normalize);
    acc_out_s{n_shuff} = f_evaluate_ens_cv(ens_out_s, firing_rate_s_norm, ens_params_s);
    fprintf('--%d', n_shuff);
end
fprintf('\nDone\n');

acc_out_full = cat(2,acc_out_s{:});
thresh = prctile(acc_out_full(:), 1);

figure; 
subplot(2,1,1); hold on;
[f_d,x_d] = ecdf(acc_out_d);
[f_s,x_s] = ecdf(acc_out_full(:));
plot(x_d, f_d, 'LineWidth', 2);
plot(x_s, f_s, 'LineWidth', 2);
line([thresh thresh], [0 1], 'Color', 'r', 'LineStyle', '--')
legend('Data', 'Shuff', '99% Thresh')
xlim([min(x_d) max(x_d)*2])
bndw = 1;
xlabel('test error');
ylabel('component fraction')
title('component test error ecdf')
subplot(2,1,2); hold on;
[f_d,x_d] = ksdensity(acc_out_d, 'Bandwidth', bndw);
[f_s,x_s] = ksdensity(acc_out_full(acc_out_full<(max(x_d)*2)), 'Bandwidth', bndw);
plot(x_d, f_d, 'LineWidth', 2);
plot(x_s, f_s, 'LineWidth', 2);
line([thresh thresh], [0 max([f_d, f_s])], 'Color', 'r', 'LineStyle', '--')
legend('Data', 'Shuff', '99% Thresh')
xlim([min(x_d) max(x_d)*2])
ylim([0 max([f_d, f_s])])
title('component test error kernel density')
xlabel('test error');
ylabel('component number')

fprintf('%.1f%% ensembles are above shuff\n', sum(x_d<thresh)/numel(x_d)*100);

%% analyze ensembles
f_plot_raster_mean(firing_rate_sm(ens_out.ord_cell,:), 1);
title('raster cell sorted');

for n_comp = 1:numel(ens_out.cells.ens_list)
    cells1 = ens_out.cells.ens_list{n_comp};
    trials1 = ens_out.trials.ens_list{n_comp};
    scores1 = ens_out.cells.ens_scores(n_comp,:);
    coeffs1 = ens_out.coeffs(cells1,n_comp);

    f_plot_ensamble_deets(firing_rate_sm, cells1, trials1, scores1, coeffs1);
    title(sprintf('%s ensamble %d; test error=%.1f', ens_params.ensamble_method, n_comp, acc_out_d(n_comp)));
end


%%

figure; imagesc(ens_out.cells.ens_scores)
ens_scores_sort = f_get_stim_trig_resp(ens_out.cells.ens_scores, stim_times, trig_window);
[ens_scores_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(ens_scores_sort, trial_types, mmn_freq, app.ops);

% freqs
for n_comp = 1:ens_params.num_comp
    figure;
    for n_tr = 1:10
        curr_tt = app.ops.context_types_all(n_tr);
        subplot(2,5,n_tr); hold on; axis tight
        tr_idx = trial_types_wctx == curr_tt;
        tr_ave = squeeze(ens_scores_sort_wctx(n_comp,:,tr_idx));
        trace_mean = mean(tr_ave,2);
        trace_sem = std(tr_ave,[],2)/sqrt(sum(tr_idx)-1);
        plot(tr_ave, 'color', [.5 .5 .5]);
        plot(trace_mean, 'r', 'LineWidth', 2);
        ylim([0 0.22])
        if n_tr == 1
            title(['comp ' num2str(n_comp)])
        end
    end
end

% ctx
ctx_trials = [18 19 20 28 29 30];
for n_comp = 1:ens_params.num_comp
    figure;
    for n_tr = 1:numel(ctx_trials)
        curr_tt = app.ops.context_types_all(ctx_trials(n_tr));
        subplot(2,3,n_tr); hold on; axis tight
        tr_idx = trial_types_wctx == curr_tt;
        tr_ave = squeeze(ens_scores_sort_wctx(n_comp,:,tr_idx));
        trace_mean = mean(tr_ave,2);
        trace_sem = std(tr_ave,[],2)/sqrt(sum(tr_idx)-1);
        plot(tr_ave, 'color', [.5 .5 .5]);
        plot(trace_mean, 'color', app.ops.context_types_all_colors2{ctx_trials(n_tr)}, 'LineWidth', 2);
        ylim([0 0.22])
        if n_tr == 1
            title(['comp ' num2str(n_comp)])
        end
    end
end

figure; hold on;
plot(ens_out.cells.ens_scores(2,:))
plot(ens_out.cells.ens_scores(15,:))

end