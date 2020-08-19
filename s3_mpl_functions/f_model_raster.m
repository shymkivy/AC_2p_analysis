function f_model_raster(data, ops)

num_cells = 50;
num_trials = 40;

%load('reliability_dd12')

%% load cell distributions

tt = [20, 30];

reliab_list = cell(numel(ops.regions_to_analyze), 1);
peak_mag_list = cell(numel(ops.regions_to_analyze), 1);
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    reliab_list1 = cell(cdata.num_dsets, numel(tt));
    peak_mag_list1 = cell(cdata.num_dsets, numel(tt));
    for n_dset = 1:cdata.num_dsets
        trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
        trial_types = cdata.trial_types_pr{n_dset};
        resp_cells = logical(cdata.peak_tuned_trials_full{n_dset});
        cell_reliability_thr = cdata.tuning_all{n_dset}.peak_tuning_full_resp.stat_pk.sig_thresh;
        cell_reliability = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_reliability;
        
        for n_tt = 1:numel(tt)
            tt1 = tt(n_tt);
            resp_cells1 = resp_cells(:,tt1);
            %
            cell_reliability1 = cell_reliability(:,tt1);
            reliab_list1{n_dset, n_tt} = cell_reliability1(resp_cells1,:);
            %
            cell_reliability_thr1 = cell_reliability_thr(:,tt1);
            cell_reliability_thr2 = cell_reliability_thr1(resp_cells1,:);
            trial_peaks1 = trial_peaks(:,trial_types == ops.context_types_all(tt1));
            trial_peaks2 = trial_peaks1(resp_cells1,:);
            peak_mag_list1{n_dset, n_tt} = trial_peaks2(trial_peaks2>cell_reliability_thr2);
        end
    end
    reliab_list{n_cond} = cat(1,reliab_list1{:});
    peak_mag_list{n_cond} = cat(1,peak_mag_list1{:});
end

peak_mag_list = cat(1,peak_mag_list{:});
% figure; histogram(peak_mag_list, 20);
% [f, x] = ksdensity(peak_mag_list);
% figure;
% plot(x,f)
% 
reliab_list = cat(1,reliab_list{:});
% x = sort(reliab_list);
% f = (1:numel(reliab_list))/numel(reliab_list);
% figure; plot(x, f)

%% make rand raster
raster = zeros(num_cells, num_trials);

% first add random activity
%cell_reliability = linspace(0.15,1,num_cells);
cell_reliability = randsample(reliab_list, num_cells);

active_trials = cell(num_cells,1);
for n_cell = 1:num_cells
    num_events = round(cell_reliability(n_cell)*num_trials);
    chosen_events = randsample(num_trials, num_events);
    raster(n_cell,chosen_events) = randsample(peak_mag_list, num_events);
    active_trials{n_cell} = chosen_events;
end

%% quick analysis of rand
[dend_order_cell, ~] = f_hcluster(raster, 'cosine', 1);
[dend_order_tr, ~] = f_hcluster(raster', 'cosine', 1);

raster_sort_cell = raster(dend_order_cell,:);
raster_sort_ctr = raster_sort_cell(:,dend_order_tr);

f1 = figure;
subplot(4,3,1);
imagesc(raster)
title('model raster');

figure(f1);
subplot(4,3,2);
imagesc(raster_sort_ctr);
title('model raster sort sort');

figure(f1);
subplot(4,3,3);
imagesc(corr(raster_sort_ctr'));
title('model raster sort sort cell corr');
axis tight equal
%% ensemble params
reliab_thresh = 0.8;
ens_size = [6 6];
cell_overlap_fraction = 0;
num_corr_trials = [8 8];

%% choose ensemble cells 
ens_cell_cell = cell(numel(ens_size),1);
cell_list = find(cell_reliability<=reliab_thresh);
for n_ens = 1:(numel(ens_size)-1)
    num_corr_cells = round(ens_size(n_ens)*cell_overlap_fraction);
    for n_ens2 = (n_ens+1):(numel(ens_size))
        chosen_cells = randsample(cell_list, num_corr_cells);
        ens_cell_cell{n_ens} = cat(1, ens_cell_cell{n_ens}, chosen_cells);
        ens_cell_cell{n_ens2} = cat(1, ens_cell_cell{n_ens2}, chosen_cells);
        cell_list(logical(sum(cell_list==chosen_cells',2))) = [];
    end
end
for n_ens = 1:numel(ens_size)
    chosen_cells = randsample(cell_list,ens_size(n_ens) - numel(ens_cell_cell{n_ens}));
    ens_cell_cell{n_ens} = cat(1, ens_cell_cell{n_ens}, chosen_cells);
    cell_list(logical(sum(cell_list==chosen_cells',2))) = [];
end

%% choose ens trials
ens_trials_cell = cell(numel(ens_size),1);
trial_list = 1:num_trials;
for n_ens = 1:numel(num_corr_trials)
    ens_trials_cell{n_ens} = randsample(trial_list, num_corr_trials(n_ens))';
    trial_list(logical(sum(trial_list == ens_trials_cell{n_ens}))) = [];
end

%% make ens patterns raster
ens_pattern = zeros(num_cells, num_trials);
for n_ens = 1:numel(ens_size)
    cells1 = ens_cell_cell{n_ens};
    trials1 = ens_trials_cell{n_ens};
    for n_cell = 1:numel(cells1)
        ens_pattern(cells1(n_cell),trials1) = 1;
    end
end

peak_mag_list_high = peak_mag_list(peak_mag_list>prctile(peak_mag_list, 30));
ens_pattern_rate = zeros(num_cells, num_trials);
ens_pattern_rate(logical(ens_pattern(:))) = randsample(peak_mag_list_high, sum(ens_pattern(:)));


%% make raster
raster_ens = raster;
num_cell_activations = sum(ens_pattern,2);
% first remove some activations from given raster
for n_cell = 1:num_cells
    if num_cell_activations(n_cell)
        active_bins = find(raster_ens(n_cell,:));
        % remove some active bins to make room for ens bins
        rem_bins = randsample(active_bins, min(numel(active_bins),num_cell_activations(n_cell)));
        raster_ens(n_cell,rem_bins) = 0;
    end
end

% clear ens frames
raster_ens(logical(ens_pattern(:))) = 0;
raster_ens = raster_ens + ens_pattern_rate;
%% add endsmbles to raster
% raster_ens = raster;
% 
% available_trials = active_trials;
% 
% for n_ens = 1:numel(ens_size)
%     chosen_trials = randsample(num_trials, num_corr_trials(n_ens));
%     ens_trials_cell{n_ens} = chosen_trials;
%     for n_cell_ind = 1:ens_size(n_ens)
%         n_cell = ens_cell_cell{n_ens}(n_cell_ind);
%         for n_tr = 1:numel(chosen_trials)
%             tr1 = chosen_trials(n_tr);
%             if ~raster_ens(n_cell,tr1)
%                 if numel(available_trials{n_cell}) > 0
%                     if numel(available_trials{n_cell})>1
%                         chosen_tr = randsample(available_trials{n_cell}, 1);
%                     else
%                         chosen_tr = available_trials{n_cell};
%                     end
%                     raster_ens(n_cell,tr1) = raster_ens(n_cell,chosen_tr);
%                     raster_ens(n_cell,chosen_tr) = 0;
%                     available_trials{n_cell}(available_trials{n_cell} == chosen_tr) = [];
%                 else
%                     raster_ens(n_cell,tr1) = randsample(peak_mag_list, 1);
%                 end
%             else
%                 available_trials{n_cell}(available_trials{n_cell} == tr1) = [];
%             end
%         end
%     end
% end

%%
ens_trials1 = zeros(num_trials,1);
for n_tr = 1:numel(ens_trials_cell)
    ens_trials1(ens_trials_cell{n_tr}) = n_tr;
end
ens_cells1 = zeros(num_cells,1);
for n_cell = 1:numel(ens_cell_cell)
    ens_cells1(ens_cell_cell{n_cell}) = n_cell;
end

%% quick analysis 2

norm1 = 1;
if norm1
    raster_clust = raster_ens-mean(raster_ens,2);
else
    raster_clust = raster_ens;
end
[dend_order_cell2, ~] = f_hcluster(raster_clust, 'cosine', 1);
[dend_order_tr2, ~] = f_hcluster(raster_clust', 'cosine', 1);

raster_ens_sort_cell = raster_ens(dend_order_cell2,:);
raster_ens_sort_ctr = raster_ens_sort_cell(:,dend_order_tr2);

ens_pattern_rate_sort_cell = ens_pattern_rate(dend_order_cell2,:);
ens_pattern_rate_sort_ctr = ens_pattern_rate_sort_cell(:,dend_order_tr2);

% % spike mags
% ens_vals = cell(numel(ens_list_cell),1);
% for n_ens1 = 1:numel(ens_list_cell)
%     ens_vals{n_ens1} = cell(ens_size(n_ens1),1);
%     for n_cell_ind1 = 1:ens_size(n_ens1)
%         ens_vals{n_ens1}{n_cell_ind1} = raster_ens(ens_list_cell{n_ens1}(n_cell_ind1),ens_trials_cell{n_ens1});
%     end
% end

ens_cells_sort1 = ens_cells1(dend_order_cell2);
ens_trials_sort1 = ens_trials1(dend_order_tr2);

% ens_list_sort_cell = cell(numel(ens_size),1);
% ens_trials_sort_cell = cell(numel(ens_size),1);

%% plot


figure(f1);
subplot(4,3,4);
imagesc(raster_ens)
title('model raster ens');
f_plot_cell_indicator(raster_ens, ens_cells1, ops);
f_plot_trial_indicator2(raster, ens_trials1, 1, ops);

figure(f1);
subplot(4,3,5);
imagesc(raster_ens_sort_ctr);
title('model raster ens sort sort');
f_plot_cell_indicator(raster_ens, ens_cells_sort1, ops);
f_plot_trial_indicator2(raster, ens_trials_sort1, 1, ops);

figure(f1);
subplot(4,3,6);
imagesc(corr(raster_ens_sort_ctr'));
title('model raster ens sort sort cell corr');
axis tight equal
f_plot_cell_indicator(corr(raster_ens_sort_ctr'), ens_cells_sort1, ops);

figure(f1);
subplot(4,3,7);
imagesc(ens_pattern_rate)
title('model ens pattern rate');
f_plot_cell_indicator(raster_ens, ens_cells1, ops);
f_plot_trial_indicator2(raster, ens_trials1, 1, ops);

figure(f1);
subplot(4,3,8);
imagesc(ens_pattern_rate_sort_ctr);
title('model ens pattern rate sort sort');
f_plot_cell_indicator(raster_ens, ens_cells_sort1, ops);
f_plot_trial_indicator2(raster, ens_trials_sort1, 1, ops);

figure(f1);
subplot(4,3,9);
imagesc(corr(ens_pattern_rate_sort_ctr'));
title('model ens pattern rate sort sort cell corr');
axis tight equal
f_plot_cell_indicator(corr(ens_pattern_rate_sort_ctr'), ens_cells_sort1, ops);

%% run real ens analysis 
% raster_ens_an = raster_ens_sort_ctr;
% ens_cells2 = ens_cells_sort1;
% ens_trials2 = ens_trials_sort1;

raster_ens_an = raster_ens;
ens_cells2 = ens_cells1;
ens_trials2 = ens_trials1;

for n_ens = 1:numel(ens_cell_cell)
    fprintf('ensemble %d cells:  ', n_ens);
    cells1 = find(ens_cells2 == n_ens);
    for ii = 1:numel(cells1)
        fprintf('%d ', cells1(ii));
    end
    fprintf('\n');
    fprintf('ensemble %d trials: ', n_ens);
    trials1 = find(ens_trials2 == n_ens);
    for ii = 1:numel(trials1)
        fprintf('%d ', trials1(ii));
    end
    fprintf('\n');
end


params.cond_name = 'fake data';
params.n_dset = 99;
params.normalize = 'norm_mean'; %'norm_full' 'norm_mean' 'none'
params.num_comps = 2;
params.plot_stuff = 0;
params.ensamble_method = 'nmf';
params.use_LR_proj = 0;
params.ensamble_extraction = 'thresh'; % 'thresh', 'clust'
params.ensamble_extraction_thresh = 'shuff';
ens_out = f_ensemble_analysis_peaks3(raster_ens_an, params, ops);

%% evaluate ens
params2.plot_stuff = 1;
eval_out = f_evaluate_ens(ens_out, ens_trials2, ens_cells2, params2);

%% plot 

raster_ens_sort_cell2 = raster_ens_an(ens_out.cell_clust.dend_order,:);
raster_ens_sort_ctr2 = raster_ens_sort_cell2(:,ens_out.trial_clust.dend_order);

clust_ident_cell_align = eval_out.clust_eval_cell.aligned_seq(ens_out.cell_clust.clust_ident+1,1);
clust_ident_trial_align = eval_out.clust_eval_tr.aligned_seq(ens_out.trial_clust.clust_ident+1,1);

figure(f1);
subplot(4,3,10);
imagesc(raster_ens_sort_ctr2)
title('Post ens analysis ground truth label');
f_plot_cell_indicator(raster_ens, ens_cells2(ens_out.cell_clust.dend_order), ops);
f_plot_trial_indicator2(raster, ens_trials2(ens_out.trial_clust.dend_order), 1, ops);

figure(f1);
subplot(4,3,11);
imagesc(raster_ens_sort_ctr2);
title('Post ens analysis discovered label');
f_plot_cell_indicator(raster_ens, clust_ident_cell_align(ens_out.cell_clust.dend_order), ops);
f_plot_trial_indicator2(raster, clust_ident_trial_align(ens_out.trial_clust.dend_order), 1, ops);

figure(f1);
subplot(4,3,12);
imagesc(corr(raster_ens_sort_ctr2'));
title('Post ens analysis cell corr');
axis tight equal
f_plot_cell_indicator(corr(raster_ens_sort_ctr'), clust_ident_cell_align(ens_out.cell_clust.dend_order), ops);

%%
disp('Done')
end

