function stat_out = f_mpl_stat_get_thresholds(trial_data, trial_types, trials_to_analyze, ops)

%%
waitbar_on = 0;
if ~exist('ops', 'var'); ops = struct; end
if ~isfield(ops, 'stat'); ops.stat = struct; end

if isfield(ops.stat, 'thresh_method')
    thresh_method = ops.stat.thresh_method;
else
    thresh_method = 'ecdf_percentile'; %'ecdf_percentile', 'zscore_around_mean', 'zscore_around_median'
end

if isfield(ops.stat, 'num_samples_drawn')
    num_samples_drawn = ops.stat.num_samples_drawn;
else
    num_samples_drawn = 1000; % default
end

if isfield(ops.stat, 'z_scores_thresh')
    z_scores_thresh = ops.stat.z_scores_thresh;
else
    z_scores_thresh = 2; % default
end

if isfield(ops.stat, 'ecdf_percentile_thresh')
    ecdf_percentile_thresh = ops.stat.ecdf_percentile_thresh;
else
    ecdf_percentile_thresh = 95; % default
end

if isfield(ops.stat, 'min_samp_size_z_normalization')
    min_samp_size_z_normalization = ops.stat.min_samp_size_z_normalization;
else
    min_samp_size_z_normalization = 1; % default
end
%%
num_dims = ndims(trial_data);
if num_dims == 2
    [num_cells, num_trials_all] = size(trial_data);
    num_bins = 1;
    trial_data = reshape(trial_data, num_cells, num_bins, num_trials_all);
else
    [num_cells, num_bins, num_trials_all] = size(trial_data);
end

trial_data2 = trial_data(:,:,logical(sum(trial_types == trials_to_analyze',2)));
trial_types2 = trial_types(logical(sum(trial_types == trials_to_analyze',2)));
[num_cells, num_bins, num_trials_all] = size(trial_data2);

if isfield(ops.stat, 'plot_examples')
    plot_examples = randsample(1:num_cells, ops.stat.plot_examples);
else
    plot_examples = 0;
end


num_trial_cat =  numel(trials_to_analyze);
num_trials_per_cat = zeros(num_trial_cat,1);
sorted_trials = cell(num_trial_cat,1);
trial_index = 1:num_trials_all;
for n_trial = 1:num_trial_cat
    sorted_trials{n_trial} = trial_index(trials_to_analyze(n_trial) == trial_types2);
    num_trials_per_cat(n_trial) = sum(trials_to_analyze(n_trial) == trial_types2);
end

if min_samp_size_z_normalization
    num_trials_per_cat = ones(num_trial_cat,1)*min(num_trials_per_cat);
end

if waitbar_on
    wb = f_waitbar_initialize([], 'Statistics...');
end
% process: get distribution of trial average trage for given sample number,
% and approximate the location of 2 or 3 std cut off (95 or 99%). can get
% the cut off empirically, or compute if from the part of the normal curve
% that is available 
sig_thresh = zeros(num_cells, num_bins, num_trial_cat);
z_factors_out = zeros(num_cells, num_bins, num_trial_cat);
means_out = zeros(num_cells, num_bins, num_trial_cat);
for n_cell = 1:num_cells
    
    if sum(n_cell == plot_examples) && ~strcmp(thresh_method, 'zscore_around_mean')
        figure; hold on;
    end    
    
    % first estimate SEM around mean
    data_mean = mean(trial_data2(n_cell,:,:),3);
    for n_trial = 1:num_trial_cat
        data_SEM_est = std(trial_data2(n_cell,:,:),[],3)/sqrt(num_trials_per_cat(n_trial)-1);
        
        [~, bin_ind] = max(data_mean+data_SEM_est);
        
        if strcmp(thresh_method, 'zscore_around_mean')
            sig_thresh(n_cell, :, n_trial) = data_mean + z_scores_thresh*data_SEM_est;
            z_factors_out(n_cell, :, n_trial) = data_SEM_est;
            means_out(n_cell, :, n_trial) = data_mean;
        else
            samp_data = zeros(num_samples_drawn, num_bins);
            subsamp_trials = zeros(num_trials_per_cat(n_trial),1);
            for n_samp = 1:num_samples_drawn
                for n_subsamp = 1:num_trials_per_cat(n_trial)
                    % because diff num of samples per cat, have to adjust probability per cat
                    cat_num = ceil(rand(1)*num_trial_cat);
                    subsamp_trials(n_subsamp) = randsample(sorted_trials{cat_num},1);
                end
                samp_data(n_samp,:) = mean(trial_data2(n_cell,:,subsamp_trials),3);
            end
            
            % compute SEM around median
            samp_data_median = median(samp_data);
            samp_data_med_sub = samp_data - samp_data_median;
            samp_data_SEM_ar_median = sqrt(sum(samp_data_med_sub.^2)./sum(logical(samp_data_med_sub)));

%             ecdf_thresh = zeros(num_bins,1);
%             for n_bin = 1:num_bins
%                 [f,x] = ecdf(samp_data(:,n_bin));
%                 f_thresh_ind = find(f >= ecdf_percentile_thresh/100);
%                 ecdf_thresh(n_bin) = x(f_thresh_ind(1));
%             end
            
            samp_data_sort = sort(samp_data);
            ecdf_thresh = samp_data_sort(round(ecdf_percentile_thresh/100*num_samples_drawn),:);
            
            if strcmp(thresh_method, 'ecdf_percentile')
                sig_thresh(n_cell, :, n_trial) = ecdf_thresh;
                z_factors_out(n_cell, :, n_trial) = (samp_data_sort(round(0.95*num_samples_drawn),:)-samp_data_median)/2;
                means_out(n_cell, :, n_trial) = samp_data_median;
            elseif strcmp(thresh_method, 'zscore_around_median')
                sig_thresh(n_cell, :, n_trial) = samp_data_median + z_scores_thresh*samp_data_SEM_ar_median;
                z_factors_out(n_cell, :, n_trial) = samp_data_SEM_ar_median;
                means_out(n_cell, :, n_trial) = samp_data_median;
            end
            
            
            if sum(n_cell == plot_examples) && ~strcmp(thresh_method, 'zscore_around_mean')
                % compute SEM around mean
                samp_data_mean = mean(samp_data);
                samp_data_mean_sub = samp_data - samp_data_mean;
                samp_data_SEM_ar_mean = sqrt(sum(samp_data_mean_sub.^2)./sum(logical(samp_data_mean_sub)));
                ecdf_thresh_95 = samp_data_sort(round(.95*num_samples_drawn),:);
                ecdf_thresh_99 = samp_data_sort(round(.99*num_samples_drawn),:);
                
                subplot(2,5,n_trial); hold on;
                plot(samp_data_sort(:,bin_ind),1/num_samples_drawn:1/num_samples_drawn:1, 'color', [0 0 1]);
                x_cdf = min(samp_data(:,bin_ind)):1/num_samples_drawn:max(samp_data(:,bin_ind));
                plot(x_cdf, cdf('Normal',x_cdf,samp_data_median(bin_ind), samp_data_SEM_ar_median(bin_ind)),'color', [1 0 1]);
                plot(x_cdf, cdf('Normal',x_cdf,samp_data_mean(bin_ind), samp_data_SEM_ar_mean(bin_ind)), 'color', [0 1 0]);
                line([samp_data_median(bin_ind) samp_data_median(bin_ind)],[0 1],'color', [1 0.5 1]);
                line([samp_data_mean(bin_ind) samp_data_mean(bin_ind)],[0 1], 'color',[0.5 1 0.5]);
                line([ecdf_thresh_95(bin_ind) ecdf_thresh_95(bin_ind)],[0 1],'color', [0.5 0.5 1],'LineStyle','--');
                line([ecdf_thresh_99(bin_ind) ecdf_thresh_99(bin_ind)],[0 1],'color', [0.5 0.5 1],'LineStyle','-.');
                line([samp_data_median(bin_ind) samp_data_median(bin_ind)]+2*samp_data_SEM_ar_median(bin_ind),[0 1],'color', [1 0.5 1],'LineStyle','--');
                line([samp_data_median(bin_ind) samp_data_median(bin_ind)]+3*samp_data_SEM_ar_median(bin_ind),[0 1],'color', [1 0.5 1],'LineStyle','-.');
                line([samp_data_mean(bin_ind) samp_data_mean(bin_ind)]+2*samp_data_SEM_ar_mean(bin_ind),[0 1],'color', [0.5 1 0.5],'LineStyle','--');
                line([samp_data_mean(bin_ind) samp_data_mean(bin_ind)]+3*samp_data_SEM_ar_mean(bin_ind),[0 1],'color', [0.5 1 0.5],'LineStyle','-.');
                line([0 0]+2*data_SEM_est(bin_ind),[0 1],'color', [1 0.5 0.5],'LineStyle','--');
                line([0 0]+3*data_SEM_est(bin_ind),[0 1],'color', [1 0.5 0.5],'LineStyle','-.');
                title(sprintf('trial %d, n=%d, bin %d', n_trial, num_trials_per_cat(n_trial),bin_ind));
            end
            
        end
    end
    if sum(n_cell == plot_examples) && ~strcmp(thresh_method, 'zscore_around_mean')
        legend('ECDF', 'CDF SEMarMed', 'CDF SEMarMean', 'median', 'mean', sprintf('%d%s ecdf thresh',95, '%'), sprintf('%d%s ecdf thresh',99.7, '%'), sprintf('SEMarMed %dz',2), sprintf('SEMarMed %dz',3), sprintf('SEMarMean %dz',2),sprintf('SEMarMean %dz',3),sprintf('est SEM %dz',2),sprintf('est SEM %dz',3),'FontSize',8);
        suptitle(sprintf('cell %d, %d samples',n_cell, num_samples_drawn))
    end
    if waitbar_on
        f_waitbar_update(wb,n_cell/num_cells);
    end
end
if waitbar_on
    f_waitbar_close(wb);
end

if num_dims == 2
    sig_thresh = reshape(sig_thresh, num_cells, []);
    z_factors_out = reshape(z_factors_out, num_cells, []);
    means_out = reshape(means_out, num_cells, []);
end

stat_out.sig_thresh = sig_thresh;
stat_out.z_factors = z_factors_out;
stat_out.means = means_out;

end