function tuning_out = f_tuning_temporal(trial_data_sort, trials_to_analyze, sig_thresh, dset_params, ops)

trial_window_t = dset_params.trial_window_t;
trial_types = dset_params.trial_types;
onset_window_frames = dset_params.onset_window_frames;
offset_window_frames = dset_params.offset_window_frames;

trial_ave_freq = f_mpl_trial_average(trial_data_sort,trial_types, trials_to_analyze, 'none');
trials_to_analyze_indx = logical(sum(trial_types == trials_to_analyze',2));

if isempty(sig_thresh)
    [sig_thresh, z_factors_out, means_out, ~] = f_mpl_stat_get_thresholds(trial_data_sort(:,:,trials_to_analyze_indx), trial_types(trials_to_analyze_indx), ops);
%             ops2 = ops; ops2.stat.thresh_method = 'zscore_around_mean';
%             [thresh_out, z_out, mean_out] = f_mpl_stat_get_thresholds(trial_data_sort_sort(:,:,1:400), trial_types(1:400), ops2);%         
end


% compute tunning 
max_onset_resp = squeeze(max(trial_ave_freq(:,onset_window_frames,:),[],2));
max_onset_thresh = squeeze(max(sig_thresh(:,onset_window_frames,:),[],2));
max_offset_resp = squeeze(max(trial_ave_freq(:,offset_window_frames,:),[],2));
max_offset_thresh = squeeze(max(sig_thresh(:,offset_window_frames,:),[],2));

% look for any points in the trace that cross the thresh
tuning_z = trial_ave_freq./sig_thresh;
tuning_z_ind = tuning_z>1;

num_cells = size(trial_ave_freq,1);
onset_tunning = zeros(num_cells, numel(trials_to_analyze));
onset_tunning_mag = zeros(num_cells, numel(trials_to_analyze));
onset_sensitivity = zeros(num_cells, numel(trials_to_analyze));
for n_cell = 1:num_cells
    for n_trial = 1:numel(trials_to_analyze)
        temp_trace = trial_ave_freq(n_cell,onset_window_frames,n_trial);
        temp_thresh = sig_thresh(n_cell,onset_window_frames,n_trial);
        onset_tunning(n_cell,n_trial) = sum(temp_trace(temp_trace == max(temp_trace))>temp_thresh(temp_trace == max(temp_trace)));
        base_mean = mean(mean(trial_ave_freq(n_cell,onset_window_frames,(1:numel(trials_to_analyze))~=n_trial),3),2);
        %onset_tunning_mag(n_cell,n_trial) = (mean(temp_trace) - mean(mean(trial_ave_freq(n_cell,onset_window_frames,:),3),2))/(mean(temp_trace) + mean(mean(trial_ave_freq(n_cell,onset_window_frames,:),3),2));
        %onset_tunning_mag(n_cell,n_trial) = (mean(temp_trace))/mean(mean(trial_ave_freq(n_cell,onset_window_frames,:),3),2);
        onset_tunning_mag(n_cell,n_trial) = (mean(temp_trace)-base_mean)/(mean(temp_trace)+base_mean);
        temp_data = squeeze(trial_data_sort(n_cell,onset_window_frames, trial_types == trials_to_analyze(n_trial)));
        onset_sensitivity(n_cell, n_trial) = sum(logical(sum(temp_data>(temp_thresh'))))/sum(trial_types == trials_to_analyze(n_trial));
    end
end
onset_tunning_mag(isnan(onset_tunning_mag)) = 0;

offset_tunning = zeros(num_cells, numel(trials_to_analyze));
offset_tunning_mag = zeros(num_cells, numel(trials_to_analyze));
offset_sensitivity = zeros(num_cells, numel(trials_to_analyze));
for n_cell = 1:num_cells
    for n_trial = 1:numel(trials_to_analyze)
        temp_trace = trial_ave_freq(n_cell,offset_window_frames,n_trial);
        temp_thresh = sig_thresh(n_cell,offset_window_frames,n_trial);
        offset_tunning(n_cell,n_trial) = sum(temp_trace(temp_trace == max(temp_trace))>temp_thresh(temp_trace == max(temp_trace)));
        base_mean = mean(mean(trial_ave_freq(n_cell,offset_window_frames,(1:numel(trials_to_analyze))~=n_trial),3),2);
        offset_tunning_mag(n_cell,n_trial) = (mean(temp_trace)-base_mean)/(mean(temp_trace)+base_mean);
        temp_data = squeeze(trial_data_sort(n_cell,offset_window_frames, trial_types == trials_to_analyze(n_trial)));
        offset_sensitivity(n_cell, n_trial) = sum(logical(sum(temp_data>(temp_thresh'))))/sum(trial_types == trials_to_analyze(n_trial));
    end
end
offset_tunning_mag(isnan(offset_tunning_mag)) = 0;

tuning_out.sig_thresh = sig_thresh;
tuning_out.z_factors_out = z_factors_out;
tuning_out.means_out = means_out;
tuning_out.onset_tunning = onset_tunning;
tuning_out.onset_tunning_mag = onset_tunning_mag;
tuning_out.onset_sensitivity = onset_sensitivity;
tuning_out.offset_tunning = offset_tunning;
tuning_out.offset_tunning_mag = offset_tunning_mag;
tuning_out.offset_sensitivity = offset_sensitivity;

%Sensitivity = (True Positive)/(True Positive + False Negative)
%Specificity = (True Negative)/(True Negative + False Positive)

if ops.stat.plot_examples
    plot_cells = sort(randsample(1:num_cells, ops.stat.plot_examples));
    for n_cell_ind = 1:numel(plot_cells)%900:910
        n_cell = plot_cells(n_cell_ind);

        figure; hold on;
        plot(max_onset_resp(n_cell,:), 'g', 'linewidth', 2);
        plot(max_offset_resp(n_cell,:), 'b', 'linewidth', 2);
        plot(max_onset_thresh(n_cell,:), '--g');
        plot(max_offset_thresh(n_cell,:), '--b');
        title(sprintf('Tuning for cell %d', n_cell));
        legend('Onset', 'Offset');
        axis tight

        figure;
        subplot(1,2,1);
        imagesc(trial_window_t,1:numel(trials_to_analyze),squeeze(trial_ave_freq(n_cell,:,:))');
        xlabel('time');
        ylabel('Freq');
        title('Trial ave')
        subplot(1,2,2);
        imagesc(trial_window_t,1:numel(trials_to_analyze),squeeze(tuning_z_ind(n_cell,:,:))');
        xlabel('time');
        title('Tuning Z thesh crossed');
        suptitle(sprintf('Tuning for cell %d', n_cell));

        figure;
        y_max = max(max(trial_data_sort(n_cell,:,:)));
        if ~y_max
            y_max = 1;
        end
        if numel(trials_to_analyze)<= 10
            n = 5;
            m = 2;
        elseif numel(trials_to_analyze) <= 15
            n = 5;
            m = 3;
        elseif numel(trials_to_analyze) <= 21
            n = 7;
            m = 3;
        elseif numel(trials_to_analyze) <=28
            n = 7;
            m = 4;
        end
        for n_trial = 1:numel(trials_to_analyze)
            subplot(m,n,n_trial)
            hold on;
            plot(trial_window_t, squeeze(trial_data_sort(n_cell,:,trial_types==trials_to_analyze(n_trial))), 'color', [0.8 0.8 0.8]);
            plot(trial_window_t, squeeze(trial_ave_freq(n_cell,:,n_trial)), 'm', 'linewidth', 2);
            plot(trial_window_t, sig_thresh(n_cell,:,n_trial), '--g');
            plot(trial_window_t, y_max*ones(1,numel(trial_window_t)).*tuning_z_ind(n_cell,:,n_trial), '--r')
            axis tight;
            ylim([0 y_max]);
            text(trial_window_t(2), y_max-y_max/5, sprintf('On resp = %d\nOff resp = %d\nOn mag = %.2f\nOff mag = %.2f\nOn sens = %.2f\nOff sens = %.2f', onset_tunning(n_cell, n_trial), offset_tunning(n_cell, n_trial), onset_tunning_mag(n_cell, n_trial), offset_tunning_mag(n_cell, n_trial), onset_sensitivity(n_cell, n_trial), offset_sensitivity(n_cell, n_trial)), 'FontSize', 8)
            title(ops.context_types_labels{n_trial});
        end
        suptitle(sprintf('Cell %d freq trial firing rate', n_cell));
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