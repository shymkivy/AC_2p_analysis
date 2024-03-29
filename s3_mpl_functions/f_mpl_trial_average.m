function [trial_ave, trial_std, num_trials, baseline] = f_mpl_trial_average(trial_data_sort,trial_types, context_types, baseline_removal_method, baseline_window_frames)
% generates trial averaged data from sortedvectorized data
%   inputs:
%       trial_data_sort(n_cell,trial_frames,trial_types)
%       trial_type(n_trials)
%       context_types
%       baseline_removal_method
%
%   outputs:
%       trial_ave(n_cell,trial_frame,trial_type_ave)
%       baseline that was removed
if ~exist('baseline_removal_method', 'var') || ~sum(strcmpi(baseline_removal_method, {'kde', 'mean', 'min', 'none'}))
    baseline_removal_method = 'kde';
end
 

[num_cells, num_frames, ~] = size(trial_data_sort);

% do trial average over all contexts
trial_ave = zeros(num_cells, num_frames, numel(context_types));
trial_std = zeros(num_cells, num_frames, numel(context_types));
num_trials = zeros(numel(context_types),1);
for n_cell = 1:num_cells
    % there are two sources of variance, the gaussian white noise, and
    % the signal variance. when whe average over trials, the white
    % noise averages out, but the signal should remain. when we average
    % over trials, we shoud not divide by std, because the white noise
    % is different for every cell. To normalize everything, we adjusted
    % that peak firing rate to 1 previously.              

    % average for each context

    for n_cntxt = 1:numel(context_types)
        ctx_trials_indx = (trial_types == context_types(n_cntxt));
        num_trials(n_cntxt) = sum(ctx_trials_indx);
        if sum(ctx_trials_indx) > 0
            trial_ave(n_cell,:,n_cntxt) = mean(trial_data_sort(n_cell,:,logical(ctx_trials_indx)),3);
            trial_std(n_cell,:,n_cntxt) = std(trial_data_sort(n_cell,:,logical(ctx_trials_indx)),[],3);
        else
            trial_ave(n_cell,:,n_cntxt) = zeros(num_frames,1);
            trial_std(n_cell,:,n_cntxt) = zeros(num_frames,1);
        end
    end


    if strcmpi(baseline_removal_method, 'kde')
        % find the baseline to subtract. Use KDE and find peak
        [f,xi] = ksdensity(reshape(trial_ave(n_cell,:,:),[],1));
        baseline = mean(xi(f == max(f)));  % added mean not sure if is ok
    elseif strcmpi(baseline_removal_method, 'mean')
        baseline = mean(mean(trial_ave(n_cell, baseline_window_frames,:)));
    elseif strcmpi(baseline_removal_method, 'min')
        baseline = min(min(trial_ave(n_cell, baseline_window_frames,:)));
    elseif strcmpi(baseline_removal_method, 'none')
        baseline = 0;
    end
    trial_ave(n_cell,:,:) = trial_ave(n_cell,:,:) - baseline;

    % plot KDE
    %                 figure
    %                 plot(xi,f);
    %                 title(sprintf('cell %d',n_cell))

    %         % plots for fun
    %         figure;
    %         for ii = 1:numel(context_types)
    %             subplot(4,7,ii)
    %             %plot(temp_data(:,temp_trial_type == context_types(ii)))
    %             hold on;
    %             plot(trial_ave(n_cell,:,ii) , 'k', 'LineWidth', 2);
    %             ylim([-0.1 0.3]);
    %             title(context_type2(ii));
    %         end
    %         suptitle(sprintf('cell %d',n_cell))
end



end