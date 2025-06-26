function f_dv_ensemble_plot(app)

ddata = app.ddata;
ens_stats = ddata.ensemble_stats{1};


if isempty(ens_stats)
    disp('Run ensemble detection and stats first')
else
    
    ensembles = ddata.ensembles{1};
    ens_params = ens_stats.ens_params;

    cdata = f_dv_get_cdata(app);
    firing_rate = cat(1,cdata.S_sm);
    active_cells = sum(firing_rate,2) ~= 0;
    firing_rate(~active_cells,:) = [];

    firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/ddata.proc_data{1}.frame_data.volume_period);
    
    acc_out_d = ensembles.acc_out_d;
    acc_out_s = ens_stats.acc_out_shuff;
    acc_out_full = cat(1,acc_out_s{:});
    thresh = ens_stats.acc_out_thresh;
   
    %%
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

    fprintf('%.1f%% (%d) ensembles are above shuff %.1f%% thresh\n', sum(acc_out_d<thresh)/numel(acc_out_d)*100, sum(acc_out_d<thresh), ens_params.shuff_thresh_percent);

    %% analyze ensembles
    f_plot_raster_mean(firing_rate_sm(ensembles.ord_cell,:), 1);
    title('raster cell sorted');

    for n_comp = 1:numel(ensembles.cells.ens_list)
        plot1 = 0;
        if app.ensonlysignificantCheckBox.Value
            if ens_stats.accepted_ensembles(n_comp)
                plot1 = 1;
            end
        else
            plot1 = 1;
        end
        cells1 = ensembles.cells.ens_list{n_comp};
        if numel(cells1) < 2
            plot1 = 0;
        end
            
        if plot1
            trials1 = ensembles.trials.ens_list{n_comp};
            scores1 = ensembles.cells.ens_scores(n_comp,:);
            coeffs1 = ensembles.coeffs(cells1,n_comp);

            f_plot_ensemble_deets(firing_rate_sm, cells1, trials1, scores1, coeffs1);
            title(sprintf('%s ensemble %d; test error=%.1f/%.1f, %.1f%%', ens_params.ensemble_method, n_comp, acc_out_d(n_comp), thresh, acc_out_d(n_comp)/thresh*100));
        end
    end


    %%
%     stim_times = ddata.stim_frame_index{n_pl};
%     mmn_freq = ddata.MMN_freq{1};
%     trial_types = ddata.trial_types{1};
% 
%     trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
%     [~, trial_frames] = f_dv_compute_window_t(trial_window, cdata.volume_period);
% 
% 
%     figure; imagesc(ens_out.cells.ens_scores)
%     ens_scores_sort = f_get_stim_trig_resp(ens_out.cells.ens_scores, stim_times, trial_frames);
%     [ens_scores_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(ens_scores_sort, trial_types, mmn_freq, app.ops);
% 
%     % freqs
%     for n_comp = 1:ens_params.num_comp
%         figure;
%         for n_tr = 1:10
%             curr_tt = app.ops.context_types_all(n_tr);
%             subplot(2,5,n_tr); hold on; axis tight
%             tr_idx = trial_types_wctx == curr_tt;
%             tr_ave = squeeze(ens_scores_sort_wctx(n_comp,:,tr_idx));
%             trace_mean = mean(tr_ave,2);
%             trace_sem = std(tr_ave,[],2)/sqrt(sum(tr_idx)-1);
%             plot(tr_ave, 'color', [.5 .5 .5]);
%             plot(trace_mean, 'r', 'LineWidth', 2);
%             ylim([0 0.22])
%             if n_tr == 1
%                 title(['comp ' num2str(n_comp)])
%             end
%         end
%     end
% 
%     % ctx
%     ctx_trials = [18 19 20 28 29 30];
%     for n_comp = 1:ens_params.num_comp
%         figure;
%         for n_tr = 1:numel(ctx_trials)
%             curr_tt = app.ops.context_types_all(ctx_trials(n_tr));
%             subplot(2,3,n_tr); hold on; axis tight
%             tr_idx = trial_types_wctx == curr_tt;
%             tr_ave = squeeze(ens_scores_sort_wctx(n_comp,:,tr_idx));
%             trace_mean = mean(tr_ave,2);
%             trace_sem = std(tr_ave,[],2)/sqrt(sum(tr_idx)-1);
%             plot(tr_ave, 'color', [.5 .5 .5]);
%             plot(trace_mean, 'color', app.ops.context_types_all_colors2{ctx_trials(n_tr)}, 'LineWidth', 2);
%             ylim([0 0.22])
%             if n_tr == 1
%                 title(['comp ' num2str(n_comp)])
%             end
%         end
%     end
% 
%     figure; hold on;
%     plot(ens_out.cells.ens_scores(2,:))
%     plot(ens_out.cells.ens_scores(15,:))
end
end