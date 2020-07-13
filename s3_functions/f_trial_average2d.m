function [trial_ave, baseline] = f_trial_average2d(trial_data_sort,trial_type,ops)
% generates trial averaged data from sortedvectorized data
%   inputs:
%       trial_data_sort(n_cell,trial_frames,trial_types)
%       trial_type(n_trials)
%       ops
%
%   outputs:
%       trial_ave(n_cell,trial_frame,trial_type_ave)
%       baseline that was removed


[d1, d2, num_frames] = size(trial_data_sort);

%[num_cells, num_frames, ~] = size(trial_data_sort);

% do trial average over all contexts
trial_ave = zeros(num_cells, num_frames, numel(ops.context_types_all));
for n_cell = 1:num_cells
    % there are two sources of variance, the gaussian white noise, and
    % the signal variance. when whe average over trials, the white
    % noise averages out, but the signal should remain. when we average
    % over trials, we shoud not divide by std, because the white noise
    % is different for every cell. To normalize everything, we adjusted
    % that peak firing rate to 1 previously.              

    % average for each context

    for n_cntxt = 1:numel(ops.context_types_all)
        ctx_trials_indx = (trial_type == ops.context_types_all(n_cntxt));
        if sum(ctx_trials_indx) > 0
            trial_ave(n_cell,:,n_cntxt) = mean(trial_data_sort(n_cell,:,ctx_trials_indx),3);
        else
            trial_ave(n_cell,:,n_cntxt) = zeros(num_frames,1);
        end
    end


    if ops.baseline_removaltrial_ave == 1
        % find the baseline to subtract. Use KDE and find peak
        [f,xi] = ksdensity(reshape(trial_ave(n_cell,:,:),[],1));
        baseline = mean(xi(f == max(f)));  % added mean not sure if is ok
    elseif ops.baseline_removaltrial_ave == 2
        baseline = mean(mean(trial_ave(n_cell,ops.win.baseline_window_frames,:)));
    elseif ops.baseline_removaltrial_ave == 3
        baseline = min(min(trial_ave(n_cell,ops.win.baseline_window_frames,:)));
    else
        baseline = 0;
    end
    trial_ave(n_cell,:,:) = trial_ave(n_cell,:,:) - baseline;

    % plot KDE
    %                 figure
    %                 plot(xi,f);
    %                 title(sprintf('cell %d',n_cell))

    %         % plots for fun
    %         figure;
    %         for ii = 1:numel(ops.context_types_all)
    %             subplot(4,7,ii)
    %             %plot(ops.time_stim_window, temp_data(:,temp_trial_type == ops.context_types_all(ii)))
    %             hold on;
    %             plot(ops.time_stim_window, trial_ave(n_cell,:,ii) , 'k', 'LineWidth', 2);
    %             ylim([-0.1 0.3]);
    %             title(context_type2(ii));
    %         end
    %         suptitle(sprintf('cell %d',n_cell))
end



end