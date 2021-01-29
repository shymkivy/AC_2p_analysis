function f_mpl_plot_cell_details(data, ops)
%ops.stat.plot_examples = 5;

if ops.stat.plot_examples
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        cdata = data(strcmpi(data.area, cond_name),:);
        num_dsets = numel(cdata.area);

        dset_cell_ind = cell(num_dsets,1);
        for n_dset = 1:num_dsets
            dset_cell_ind{n_dset} = [ones(cdata.num_cells(n_dset),1)*n_dset,...
                                (1:cdata.num_cells(n_dset))'];
        end
        dset_cell_ind = cat(1,dset_cell_ind{:});

        pk_tuned_trials = cat(1,cdata.peak_tuned_trials_combined_ctx{:});
        pk_tuned_cell_ind = find(logical(sum(pk_tuned_trials,2)));

        plot_cells = sort(randsample(pk_tuned_cell_ind, ops.stat.plot_examples));
        for n_cell_ind = 1:ops.stat.plot_examples%21:30%900:910
            n_dset = dset_cell_ind(plot_cells(n_cell_ind),1);
            n_cell = dset_cell_ind(plot_cells(n_cell_ind),2);
            trial_data_sort = cdata.trial_data_sort_sm_wctx{n_dset}(n_cell,:,:);
            trial_types = cdata.trial_types_wctx{n_dset};
            ctx_mmn = cdata.ctx_mmn{n_dset};
            tuning_all = cdata.tuning_all{n_dset};
            trial_window_t = cdata.trial_window{n_dset}.trial_window_t;

            [pk_resp_mag_full,pk_resp_trial_ind_full] = max(tuning_all.peak_tuning_full_resp.fr_peak_mag_ave_z(n_cell,ctx_mmn));
            pk_resp_rel_full = tuning_all.peak_tuning_full_resp.fr_peak_reliability(n_cell,ctx_mmn);
            [pk_resp_mag_on,pk_resp_trial_ind_on] = max(tuning_all.peak_tuning_onset.fr_peak_mag_ave_z(n_cell,ctx_mmn));
            pk_resp_rel_on = tuning_all.peak_tuning_onset.fr_peak_reliability(n_cell,ctx_mmn);
            [pk_resp_mag_off,pk_resp_trial_ind_off] = max(tuning_all.peak_tuning_offset.fr_peak_mag_ave_z(n_cell,ctx_mmn));
            pk_resp_rel_off = tuning_all.peak_tuning_offset.fr_peak_reliability(n_cell,ctx_mmn);

            max_onset_resp = squeeze(max(tuning_all.trace_tuning.trial_ave(n_cell,cdata.trial_window{n_dset}.onset_window_frames,:),[],2));
            max_onset_thresh = squeeze(max(tuning_all.trace_tuning.stat_trace.sig_thresh(n_cell,cdata.trial_window{n_dset}.onset_window_frames,:),[],2));
            max_offset_resp = squeeze(max(tuning_all.trace_tuning.trial_ave(n_cell,cdata.trial_window{n_dset}.offset_window_frames,:),[],2));
            max_offset_thresh = squeeze(max(tuning_all.trace_tuning.stat_trace.sig_thresh(n_cell,cdata.trial_window{n_dset}.offset_window_frames,:),[],2));

            tuning_ind = tuning_all.trace_tuning.trial_ave>tuning_all.trace_tuning.stat_trace.sig_thresh;

            figure; 
            subplot(3,1,1); hold on; axis tight;
            plot(tuning_all.peak_tuning_full_resp.fr_peak_mag_ave(n_cell,:));
            plot(tuning_all.peak_tuning_full_resp.stat_pk.sig_thresh(n_cell,:), '--');
            for n_tr = 1:numel(ops.context_types_all)
                num_tr = sum(trial_types == ops.context_types_all(n_tr));
                scatter(n_tr+(1:num_tr)/num_tr/2-0.25,tuning_all.peak_tuning_full_resp.fr_peak_mag(n_cell, trial_types == ops.context_types_all(n_tr)));
            end
            legend('trial ave', 'sig thresh');
            title(sprintf('Cell %d, full resp peaks, tr=%d, zmag=%.2f, sens=%.2f', n_cell, ctx_mmn(pk_resp_trial_ind_full), pk_resp_mag_full, pk_resp_rel_full(pk_resp_trial_ind_full)));
            subplot(3,1,2); hold on; axis tight;
            plot(tuning_all.peak_tuning_onset.fr_peak_mag_ave(n_cell,:));
            plot(tuning_all.peak_tuning_onset.stat_pk.sig_thresh(n_cell,:), '--');
            for n_tr = 1:numel(ops.context_types_all)
                num_tr = sum(trial_types == ops.context_types_all(n_tr));
                scatter(n_tr+(1:num_tr)/num_tr/2-0.25, tuning_all.peak_tuning_onset.fr_peak_mag(n_cell, trial_types == ops.context_types_all(n_tr)));
            end
            title(sprintf('Onset resp peaks, tr=%d, zmag=%.2f, sens=%.2f', ctx_mmn(pk_resp_trial_ind_on), pk_resp_mag_on, pk_resp_rel_on(pk_resp_trial_ind_on)));
            subplot(3,1,3); hold on; axis tight;
            plot(tuning_all.peak_tuning_offset.fr_peak_mag_ave(n_cell,:));
            plot(tuning_all.peak_tuning_offset.stat_pk.sig_thresh(n_cell,:), '--');
            for n_tr = 1:numel(ops.context_types_all)
                num_tr = sum(trial_types == ops.context_types_all(n_tr));
                scatter(n_tr+(1:num_tr)/num_tr/2-0.25, tuning_all.peak_tuning_offset.fr_peak_mag(n_cell, trial_types == ops.context_types_all(n_tr)));
            end
            title(sprintf('Offset resp peaks, tr=%d, zmag=%.2f, sens=%.2f', ctx_mmn(pk_resp_trial_ind_off), pk_resp_mag_off, pk_resp_rel_off(pk_resp_trial_ind_off)));


            figure; 
            subplot(2,2,1:2);hold on;
            plot(max_onset_resp, 'g', 'linewidth', 2);
            plot(max_offset_resp, 'b', 'linewidth', 2);
            plot(max_onset_thresh, '--g');
            plot(max_offset_thresh, '--b');
            title(sprintf('Tuning for cell %d', n_cell));
            legend('Onset', 'Offset');
            axis tight
            subplot(2,2,3);
            imagesc(cdata.trial_window{n_dset}.trial_window_t,1:numel(ops.context_types_all),squeeze(tuning_all.trace_tuning.trial_ave(n_cell,:,:))');
            xlabel('time');
            ylabel('Freq');
            title('Trial ave')
            subplot(2,2,4);
            imagesc(cdata.trial_window{n_dset}.trial_window_t,1:numel(ops.context_types_all),squeeze(tuning_ind(n_cell,:,:))');
            xlabel('time');
            title('Tuning Z thesh crossed');

            figure;
            y_max = max(max(trial_data_sort(:,:,:)));
            if ~y_max; y_max = 1; end
            if numel(ops.context_types_all)<= 10; n = 5; m = 2;
            elseif numel(ops.context_types_all) <= 15; n = 5; m = 3;
            elseif numel(ops.context_types_all) <= 21; n = 7; m = 3;
            elseif numel(ops.context_types_all) <=28; n = 7; m = 4;
            elseif numel(ops.context_types_all) <=30; n = 5; m = 6;
            end
            for n_trial = 1:numel(ops.context_types_all)
                subplot(m,n,n_trial)
                hold on;
                plot(trial_window_t, squeeze(trial_data_sort(:,:,trial_types==ops.context_types_all(n_trial))), 'color', [0.8 0.8 0.8]);
                plot(trial_window_t, squeeze(tuning_all.trace_tuning.trial_ave(n_cell,:,n_trial)), 'm', 'linewidth', 2);
                plot(trial_window_t, tuning_all.trace_tuning.stat_trace.sig_thresh(n_cell,:,n_trial), '--g');
                plot(trial_window_t, y_max*ones(1,numel(trial_window_t)).*tuning_ind(n_cell,:,n_trial), '--r')
                axis tight;
                ylim([0 y_max]);
                text(trial_window_t(2), y_max-y_max/5, sprintf('On resp = %d\nOff resp = %d\nOn mag = %.2f\nOff mag = %.2f\nOn sens = %.2f\nOff sens = %.2f', tuning_all.trace_tuning.onset_tuned_trials(n_cell, n_trial), tuning_all.trace_tuning.offset_tuned_trials(n_cell, n_trial), tuning_all.trace_tuning.onset_tunning_metric(n_cell, n_trial), tuning_all.trace_tuning.offset_tunning_metric(n_cell, n_trial), tuning_all.trace_tuning.onset_sensitivity(n_cell, n_trial), tuning_all.trace_tuning.offset_sensitivity(n_cell, n_trial)), 'FontSize', 8)
                if n_trial == 1
                    title(sprintf('%s, %s, Cell %d, %s', cond_name, ops.file_names.(cond_name){n_dset}, n_cell,ops.context_types_labels{n_trial}), 'Interpreter', 'none');
                elseif n_trial == 20
                    title([ops.context_types_labels{n_trial} ' freq ' num2str(cdata.MMN_freq{n_dset}(2))]);
                elseif n_trial == 30
                    title([ops.context_types_labels{n_trial} ' freq ' num2str(cdata.MMN_freq{n_dset}(1))]);
                else
                    title(ops.context_types_labels{n_trial});
                end
            end


            ctx_mmn2 = reshape(ctx_mmn,3,2);
            %n_cell = 176;
            y_max = max(max(trial_data_sort(:,:,:)));
            figure;
            for n_tr_ind = 1:3
                if pk_resp_trial_ind_full < 4
                    n_tr = ctx_mmn2(n_tr_ind,1);
                else
                    n_tr = ctx_mmn2(n_tr_ind,2);
                end
                subplot(1,3,n_tr_ind); hold on;
                plot(trial_window_t, squeeze(trial_data_sort(:,:,trial_types==ops.context_types_all(n_tr))), 'color', [0.8 0.8 0.8]);
                plot(trial_window_t, squeeze(tuning_all.trace_tuning.trial_ave(n_cell,:,n_tr)), 'color', ops.context_colors{n_tr_ind}, 'linewidth', 2);
                plot(trial_window_t, tuning_all.trace_tuning.stat_trace.sig_thresh(n_cell,:,n_tr), '--g');
                axis tight;
                ylim([0 y_max]);

                if n_tr_ind == 1
                    title(sprintf('%s, %s, Cell %d\nOn mag=%.2f, Off mag=%.2f\nOn sens=%.2f, Off sens=%.2f', cond_name, ops.file_names.(cond_name){n_dset}, n_cell, tuning_all.trace_tuning.onset_tunning_metric(n_cell, n_tr), tuning_all.trace_tuning.offset_tunning_metric(n_cell, n_tr), tuning_all.trace_tuning.onset_sensitivity(n_cell, n_tr), tuning_all.trace_tuning.offset_sensitivity(n_cell, n_tr)), 'Interpreter', 'none');
                else
                    title(sprintf('On mag=%.2f, Off mag=%.2f\nOn sens=%.2f, Off sens=%.2f', tuning_all.trace_tuning.onset_tunning_metric(n_cell, n_tr), tuning_all.trace_tuning.offset_tunning_metric(n_cell, n_tr), tuning_all.trace_tuning.onset_sensitivity(n_cell, n_tr), tuning_all.trace_tuning.offset_sensitivity(n_cell, n_tr)));
                end
            end


            %legend('raw trials', 'trial ave', 'ecdf thresh', 'zscore thresh', 'FontSize', 8)

    %             figure;
    %             y_max = max(max(trial_raw_data_sort(n_cell,:,1:400)));
    %             y_min = min(min(trial_raw_data_sort(n_cell,:,1:400)));
    %             for n_trial = 1:10
    %                 subplot(2,5,n_trial)
    %                 hold on;
    %                 plot(trial_window_t, squeeze(trial_raw_data_sort(n_cell,:,trial_types==n_trial)), 'color', [0.8 0.8 0.8]);
    %                 plot(trial_window_t, squeeze(trial_ave_freq(n_cell,:,n_trial))*10, 'm', 'linewidth', 2)
    %                 plot(trial_window_t, squeeze(trial_raw_ave_freq(n_cell,:,n_trial)), 'g', 'linewidth', 2)
    %                 axis tight;
    %                 ylim([y_min y_max]);
    %                 title(sprintf('Freq %d', n_trial));
    %             end
    %             suptitle(sprintf('Cell %d freq trial raw vs firing rate', n_cell));
        end
    end

end


end