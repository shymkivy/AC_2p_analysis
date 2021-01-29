function tuning_out = f_get_tuning_dset(trial_data_sort, trials_to_analyze, dset_params, ops)

trial_window_t = dset_params.trial_window{1}.trial_window_t;
trial_types = dset_params.trial_types_wctx{1};

%% extract peaks and lat info

peak_tuning_full_resp = f_get_peak_tuning(trial_data_sort, trial_types, trials_to_analyze, trial_window_t, ops.resp_window_time, ops);
peak_tuning_onset = f_get_peak_tuning(trial_data_sort, trial_types, trials_to_analyze, trial_window_t, ops.onset_window, ops);
peak_tuning_offset = f_get_peak_tuning(trial_data_sort, trial_types, trials_to_analyze, trial_window_t, ops.offset_window, ops);

% figure;histogram(fr_peak_mag(1,fr_peak_mag(1,:)>0))
% figure;histogram(fr_peak_latency(:,fr_peak_latency(:,:)>1))

%% traces

trace_tuning = f_get_trace_tuning(trial_data_sort,trial_types, trials_to_analyze,dset_params, ops);

%% save
tuning_out.peak_tuning_full_resp = peak_tuning_full_resp;
tuning_out.peak_tuning_onset = peak_tuning_onset;
tuning_out.peak_tuning_offset = peak_tuning_offset;
tuning_out.trace_tuning = trace_tuning;

% 
% tuning_out.peak_tuned_trials_onset = peak_tuning_onset.fr_peak_mag_tuned_trials;
% tuning_out.peak_tuned_trials_onset_ctx = peak_tuning_onset.fr_peak_mag_tuned_trials(:,dset_params.ctx_mmn);
% tuning_out.peak_tuned_trials_offset = peak_tuning_offset.fr_peak_mag_tuned_trials;
% tuning_out.peak_tuned_trials_offset_ctx = peak_tuning_offset.fr_peak_mag_tuned_trials(:,dset_params.ctx_mmn);
% tuning_out.peak_tuned_trials_full = peak_tuning_full_resp.fr_peak_mag_tuned_trials;
% tuning_out.peak_tuned_trials_full_ctx = peak_tuning_full_resp.fr_peak_mag_tuned_trials(:,dset_params.ctx_mmn);
% tuning_out.peak_tuned_trials_combined = fr_pk_tuned_trials_combined;
% tuning_out.peak_tuned_trials_combined_ctx = fr_pk_tuned_trials_combined_ctx;
% 

%%
%if isempty(sig_thresh)
%             ops2 = ops; ops2.stat.thresh_method = 'zscore_around_mean';
%             [thresh_out, z_out, mean_out] = f_mpl_stat_get_thresholds(trial_data_sort_sort(:,:,1:400), trial_types(1:400), ops2);%         
%end

%Sensitivity = (True Positive)/(True Positive + False Negative)
%Specificity = (True Negative)/(True Negative + False Positive)




if ops.stat.plot_examples
    ctx_mmn = dset_params.ctx_mmn{1};
    
    num_trials_per_cat = zeros(numel(trials_to_analyze),1);
    for n_tr = 1:numel(trials_to_analyze)
        num_trials_per_cat(n_tr) = sum(trial_types == trials_to_analyze(n_tr));
    end
    
    max_onset_resp = squeeze(max(trace_tuning.trial_ave(:,dset_params.trial_window{1}.onset_window_frames,:),[],2));
    max_onset_thresh = squeeze(max(trace_tuning.stat_trace.sig_thresh(:,dset_params.trial_window{1}.onset_window_frames,:),[],2));
    max_offset_resp = squeeze(max(trace_tuning.trial_ave(:,dset_params.trial_window{1}.offset_window_frames,:),[],2));
    max_offset_thresh = squeeze(max(trace_tuning.stat_trace.sig_thresh(:,dset_params.trial_window{1}.offset_window_frames,:),[],2));

    % look for any points in the trace that cross the thresh
    tuning_ind = trace_tuning.trial_ave>trace_tuning.stat_trace.sig_thresh;
    %tuned_cells_ind = find(logical(peak_tuning_offset.fr_peak_mag_tuned_cells+peak_tuning_onset.fr_peak_mag_tuned_cells));

    tuned_cell_ind = find(logical(sum(peak_tuning_full_resp.fr_peak_mag_tuned_trials,2)));

    [pk_resp_mag_full,pk_resp_trial_ind_full] = max(peak_tuning_full_resp.fr_peak_mag_ave_z(:,ctx_mmn),[],2);
    pk_resp_rel_full = peak_tuning_full_resp.fr_peak_reliability(:,ctx_mmn);
    [pk_resp_mag_on,pk_resp_trial_ind_on] = max(peak_tuning_onset.fr_peak_mag_ave_z(:,ctx_mmn),[],2);
    pk_resp_rel_on = peak_tuning_onset.fr_peak_reliability(:,ctx_mmn);
    [pk_resp_mag_off,pk_resp_trial_ind_off] = max(peak_tuning_offset.fr_peak_mag_ave_z(:,ctx_mmn),[],2);
    pk_resp_rel_off = peak_tuning_offset.fr_peak_reliability(:,ctx_mmn);
    
    
    
    plot_cells = sort(randsample(tuned_cell_ind, ops.stat.plot_examples));
    for n_cell_ind = 1:ops.stat.plot_examples%21:30%900:910
        n_cell = plot_cells(n_cell_ind);%tuned_cell_ind(n_cell_ind);
        
        figure; 
        subplot(3,1,1); hold on; axis tight;
        plot(peak_tuning_full_resp.fr_peak_mag_ave(n_cell,:));
        plot(peak_tuning_full_resp.stat_pk.sig_thresh(n_cell,:), '--');
        for n_tr = 1:numel(trials_to_analyze)
            num_tr = sum(trial_types == trials_to_analyze(n_tr));
            scatter(n_tr+(1:num_tr)/num_tr/2-0.25, peak_tuning_full_resp.fr_peak_mag(n_cell, trial_types == trials_to_analyze(n_tr)));
        end
        legend('trial ave', 'sig thresh')
        title(sprintf('Cell %d, full resp peaks, tr=%d, zmag=%.2f, sens=%.2f', n_cell, ctx_mmn(pk_resp_trial_ind_full(n_cell)), pk_resp_mag_full(n_cell), pk_resp_rel_full(n_cell,pk_resp_trial_ind_full(n_cell))));
        subplot(3,1,2); hold on; axis tight;
        plot(peak_tuning_onset.fr_peak_mag_ave(n_cell,:));
        plot(peak_tuning_onset.stat_pk.sig_thresh(n_cell,:), '--');
        for n_tr = 1:numel(trials_to_analyze)
            num_tr = sum(trial_types == trials_to_analyze(n_tr));
            scatter(n_tr+(1:num_tr)/num_tr/2-0.25, peak_tuning_onset.fr_peak_mag(n_cell, trial_types == trials_to_analyze(n_tr)));
        end
        title(sprintf('Onset resp peaks, tr=%d, zmag=%.2f, sens=%.2f', ctx_mmn(pk_resp_trial_ind_on(n_cell)), pk_resp_mag_on(n_cell), pk_resp_rel_on(n_cell,pk_resp_trial_ind_on(n_cell))));
        subplot(3,1,3); hold on; axis tight;
        plot(peak_tuning_offset.fr_peak_mag_ave(n_cell,:));
        plot(peak_tuning_offset.stat_pk.sig_thresh(n_cell,:), '--');
        for n_tr = 1:numel(trials_to_analyze)
            num_tr = sum(trial_types == trials_to_analyze(n_tr));
            scatter(n_tr+(1:num_tr)/num_tr/2-0.25, peak_tuning_offset.fr_peak_mag(n_cell, trial_types == trials_to_analyze(n_tr)));
        end
        title(sprintf('Offset resp peaks, tr=%d, zmag=%.2f, sens=%.2f', ctx_mmn(pk_resp_trial_ind_off(n_cell)), pk_resp_mag_off(n_cell), pk_resp_rel_off(n_cell,pk_resp_trial_ind_off(n_cell))));
        
        
        figure; 
        subplot(2,2,1:2);hold on;
        plot(max_onset_resp(n_cell,:), 'g', 'linewidth', 2);
        plot(max_offset_resp(n_cell,:), 'b', 'linewidth', 2);
        plot(max_onset_thresh(n_cell,:), '--g');
        plot(max_offset_thresh(n_cell,:), '--b');
        title(sprintf('Tuning for cell %d', n_cell));
        legend('Onset', 'Offset');
        axis tight
        subplot(2,2,3);
        imagesc(trial_window_t,1:numel(trials_to_analyze),squeeze(trace_tuning.trial_ave(n_cell,:,:))');
        xlabel('time');
        ylabel('Freq');
        title('Trial ave')
        subplot(2,2,4);
        imagesc(trial_window_t,1:numel(trials_to_analyze),squeeze(tuning_ind(n_cell,:,:))');
        xlabel('time');
        title('Tuning Z thesh crossed');

        figure;
        y_max = max(max(trial_data_sort(n_cell,:,:)));
        if ~y_max; y_max = 1; end
        if numel(trials_to_analyze)<= 10; n = 5; m = 2;
        elseif numel(trials_to_analyze) <= 15; n = 5; m = 3;
        elseif numel(trials_to_analyze) <= 21; n = 7; m = 3;
        elseif numel(trials_to_analyze) <=28; n = 7; m = 4;
        elseif numel(trials_to_analyze) <=30; n = 5; m = 6;
        end
        for n_trial = 1:numel(trials_to_analyze)
            subplot(m,n,n_trial)
            hold on;
            if num_trials_per_cat(n_trial)
                plot(trial_window_t, squeeze(trial_data_sort(n_cell,:,trial_types==trials_to_analyze(n_trial))), 'color', [0.8 0.8 0.8]);
                plot(trial_window_t, squeeze(trace_tuning.trial_ave(n_cell,:,n_trial)), 'm', 'linewidth', 2);
                plot(trial_window_t, trace_tuning.stat_trace.sig_thresh(n_cell,:,n_trial), '--g');
                plot(trial_window_t, y_max*ones(1,numel(trial_window_t)).*tuning_ind(n_cell,:,n_trial), '--r')
                axis tight;
                ylim([0 y_max]);
                text(trial_window_t(2), y_max-y_max/5, sprintf('On resp = %d\nOff resp = %d\nOn mag = %.2f\nOff mag = %.2f\nOn sens = %.2f\nOff sens = %.2f', trace_tuning.onset_tuned_trials(n_cell, n_trial), trace_tuning.offset_tuned_trials(n_cell, n_trial), trace_tuning.onset_tunning_metric(n_cell, n_trial), trace_tuning.offset_tunning_metric(n_cell, n_trial), trace_tuning.onset_sensitivity(n_cell, n_trial), trace_tuning.offset_sensitivity(n_cell, n_trial)), 'FontSize', 8)
                if n_trial == 1
                    title(sprintf('%s, Cell %d, %s', dset_params.experiment{1}, n_cell,ops.context_types_labels{n_trial}), 'Interpreter', 'none');
                elseif n_trial == 20
                    title([ops.context_types_labels{n_trial} ' freq ' num2str(dset_params.MMN_freq{1}(2))]);
                elseif n_trial == 30
                    title([ops.context_types_labels{n_trial} ' freq ' num2str(dset_params.MMN_freq{1}(1))]);
                else
                    title(ops.context_types_labels{n_trial});
                end
            end
        end
        
        
        ctx_mmn2 = reshape(ctx_mmn,3,2);
        %n_cell = 176;
        y_max = max(max(trial_data_sort(n_cell,:,:)));
        figure;
        for n_tr_ind = 1:3
            if pk_resp_trial_ind_full(n_cell) < 4
                n_tr = ctx_mmn2(n_tr_ind,1);
            else
                n_tr = ctx_mmn2(n_tr_ind,2);
            end
            subplot(1,3,n_tr_ind); hold on;
            plot(trial_window_t, squeeze(trial_data_sort(n_cell,:,trial_types==trials_to_analyze(n_tr))), 'color', [0.8 0.8 0.8]);
            plot(trial_window_t, squeeze(trace_tuning.trial_ave(n_cell,:,n_tr)), 'color', ops.context_colors{n_tr_ind}, 'linewidth', 2);
            plot(trial_window_t, trace_tuning.stat_trace.sig_thresh(n_cell,:,n_tr), '--g');
            axis tight;
            ylim([0 y_max]);
            
            if n_tr_ind == 1
                title(sprintf('%s, Cell %d\nOn mag=%.2f, Off mag=%.2f\nOn sens=%.2f, Off sens=%.2f', dset_params.experiment{1}, n_cell, trace_tuning.onset_tunning_metric(n_cell, n_tr), trace_tuning.offset_tunning_metric(n_cell, n_tr), trace_tuning.onset_sensitivity(n_cell, n_tr), trace_tuning.offset_sensitivity(n_cell, n_tr)), 'Interpreter', 'none');
            else
                title(sprintf('On mag=%.2f, Off mag=%.2f\nOn sens=%.2f, Off sens=%.2f', trace_tuning.onset_tunning_metric(n_cell, n_tr), trace_tuning.offset_tunning_metric(n_cell, n_tr), trace_tuning.onset_sensitivity(n_cell, n_tr), trace_tuning.offset_sensitivity(n_cell, n_tr)));
            end
        end
    end
end



end


