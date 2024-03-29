function ens_tuning = f_dv_ensemble_tuning_stats(app, params)

ens_out = params.ensembles.ens_out;
num_ens = ens_out.num_comps;

if isfield(params, 'ensemble_stats')
    accepted_ens = params.ensemble_stats.accepted_ensembles;
else
    accepted_ens = ones(num_ens, 1);
end

ens_scores = ens_out.scores(accepted_ens,:);

%% params
peak_stats = 'shuff_pool'; % 'shuff_pool', 'shuff_locwise', 'z_thresh'
peak_bin_time = .250; % sec

num_samp = 2000;
%stat_resp_window = [.05 1];
%stat_resp_window = [-2 3];

z_thresh = params.z_thresh_new;
peak_prcntle = normcdf(params.z_thresh_new)*100;

%%
n_pl = params.n_pl;
ddata = params.ddata;

num_ens = size(ens_scores,1);
stim_times = cat(2,ddata.stim_frame_index{params.planes});
stim_times2 = round(mean(stim_times,2));
trial_types = ddata.trial_types{1};
MMN_freq = ddata.MMN_freq{1};
fr = 1000/double(ddata.proc_data{1}.frame_data.volume_period);
peak_bin_size = ceil(peak_bin_time*fr);

stat_window = f_str_to_array(app.stats_BaserespwinEditField.Value);
[stat_window_t, stat_frames] = f_dv_compute_window_t(stat_window, ddata.proc_data{1}.frame_data.volume_period_ave);


%%
win1 = stat_frames;

trial_data_sort = f_get_stim_trig_resp(ens_scores, stim_times2, win1);

if ~isempty(MMN_freq)
    [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, MMN_freq, app.ops);
else
    trial_data_sort_wctx = trial_data_sort;
    trial_types_wctx = trial_types;
end
% 
%% choose population for shuffle
if strcmpi(params.stat_source, 'All')
    pop_stim_times = stim_times;
    trial_data_sort_stat = trial_data_sort;
elseif strcmpi(params.stat_source, 'Freqs')
    stim_idx = logical(sum(trial_types == 1:app.ops.stim.num_freqs,2));
    pop_stim_times = stim_times(stim_idx);
    trial_data_sort_stat = trial_data_sort(:,:,stim_idx);
end

num_tt = numel(app.ops.context_types_all);
num_trials = numel(pop_stim_times);
num_trial_per_stim = num_trials/app.ops.stim.num_freqs;
num_t = sum(stat_frames);

trial_data_sort_stat_mean = squeeze(mean(trial_data_sort_stat,2));
%% convert to z scores

stat_trials_mean = mean(trial_data_sort_stat,3);
stat_trials_sem = std(trial_data_sort_stat,[],3)/sqrt(num_trial_per_stim-1);

%% get peak resp

peak_vals = zeros(num_ens, num_tt);
peak_locs = zeros(num_ens, num_tt);
for n_tt = 1:num_tt
    trial_data_sort2 = trial_data_sort_wctx(:,:, trial_types_wctx==app.ops.context_types_all(n_tt));
    [peak_vals(:,n_tt), peak_locs(:,n_tt)] = f_get_trial_peak(trial_data_sort2, peak_bin_size);
end

if strcmpi(peak_stats, 'shuff_pool') || strcmpi(peak_stats, 'shuff_locwise')
    % compute stats for peaks
    samp_peak_vals = zeros(num_ens, num_samp);
    samp_peak_locs = zeros(num_ens, num_samp);
    wb = f_waitbar_initialize(app, 'Stats: sampling...');
    for n_cell = 1:num_ens
        f_waitbar_update(wb, n_cell/num_ens, 'Stats: sampling...');
        for n_samp = 1:num_samp
            samp_idx = randsample(num_trials, num_trial_per_stim, 1);
            samp_trial_data_sort = trial_data_sort_wctx(n_cell, :, samp_idx);
            [samp_peak_vals(n_cell,n_samp), samp_peak_locs(n_cell,n_samp)] = f_get_trial_peak(samp_trial_data_sort, peak_bin_size);
        end
        %fprintf('%d-', n_cell)
    end
    f_waitbar_close(wb);
    %fprintf('\n')
end

% figure; histogram(peak_locs(:))
% figure; histogram(samp_peak_locs(n_cell,:))
% figure; histogram(samp_peak_vals(n_cell,:))

%% get resp cell

if strcmpi(peak_stats, 'shuff_pool')
    peak_resp_thresh = repmat(prctile(samp_peak_vals', peak_prcntle)', [1 num_t]);
    peak_resp_cells = peak_vals>peak_resp_thresh(:,1);
elseif strcmpi(peak_stats, 'shuff_locwise')
    peak_resp_thresh = zeros(num_ens, num_t);
    peak_resp_cells = zeros(num_ens, num_tt);
    for n_cell = 1:num_ens
        for n_loc = 1:num_t
            temp_th = prctile(samp_peak_vals(n_cell,samp_peak_locs(n_cell,:) == n_loc), peak_prcntle);
            if isnan(temp_th)
                peak_resp_thresh(n_cell, n_loc) = 0;
            else
                peak_resp_thresh(n_cell, n_loc) = temp_th;
            end
        end
        for n_tt = 1:num_tt
            peak_resp_cells(n_cell, n_tt) = peak_vals(n_cell, n_tt) > peak_resp_thresh(n_cell,peak_locs(n_cell,n_tt));
        end
    end
else
    peak_resp_thresh = zeros(num_ens, num_t);
    peak_resp_cells = zeros(num_ens, num_tt);
    for n_cell = 1:num_ens
        trial_data1 = squeeze(trial_data_sort_stat(n_cell,:,:));
        peak_resp_thresh(n_cell,:) = mean(trial_data1,2) + z_thresh*std(trial_data1,[],2)/sqrt(num_trial_per_stim-1);
        for n_tt = 1:num_tt
            peak_resp_cells(n_cell, n_tt) = peak_vals(n_cell, n_tt) > peak_resp_thresh(n_cell,peak_locs(n_cell,n_tt));
        end
    end
end

%%

volt_dat = params.ddata.proc_data{1}.volt_data_binned{n_pl};

loco1 = volt_dat(:,3);
firing_rate = ens_scores;

cell_corr = corr(loco1, firing_rate');

num_shuff = 50;
samp_data = zeros(num_shuff, num_ens);

for n_shuff = 1:num_shuff
    firing_rate_s = f_shuffle_data(firing_rate);
    samp_data(n_shuff,:) = corr(loco1, firing_rate_s');
end
loco_thresh = prctile(samp_data(:), 99);

% [f_d,x_d] = ecdf(cell_corr);
% [f_s,x_s] = ecdf(samp_data(:));
% figure; hold on;
% plot(x_d, f_d, 'LineWidth', 2);
% plot(x_s, f_s, 'LineWidth', 2);
% l1 = line([loco_thresh, loco_thresh], [0 1]);
% l1.LineStyle = '--';
% l1.Color = 'r';
% legend('loco data', 'shuff', '99% thresh');

z_factor = std(samp_data(:));

loco_cell = cell_corr>loco_thresh;
loco_corr = cell_corr;
loco_z = cell_corr/z_factor;

%%
ens_tuning.stat_trials_mean = stat_trials_mean;
ens_tuning.stat_trials_sem = stat_trials_sem;
ens_tuning.stat_trials_mean_mean = mean(stat_trials_mean,2);
ens_tuning.stat_trials_mean_sem = mean(stat_trials_sem,2);
ens_tuning.peak_resp_cells = peak_resp_cells;
ens_tuning.peak_resp_thresh = peak_resp_thresh;
ens_tuning.peak_val = peak_vals;
ens_tuning.peak_t_all = stat_window_t(peak_locs);
ens_tuning.stat_window_t = stat_window_t; % stat_window_t
ens_tuning.z_thresh = z_thresh;
ens_tuning.peak_stats = peak_stats;
ens_tuning.peak_bin_size = peak_bin_size;
ens_tuning.peak_prcntle = peak_prcntle;
ens_tuning.num_samp = num_samp;
ens_tuning.num_cells = num_ens;
ens_tuning.accepted_cells = params.cdata.accepted_cells;
ens_tuning.loco_cell = loco_cell;
ens_tuning.loco_corr = loco_corr;
ens_tuning.loco_z = loco_z;

end