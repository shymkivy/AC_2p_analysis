function f_model_raster(data, ops)

num_cells = 50;
num_trials = 40;

%load('reliability_dd12')

%% load cell traces

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

%%
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

figure;
imagesc(raster)
title('model raster');

%% analysis
[dend_order_cell, ~] = f_hcluster(raster, 'cosine', 1);
[dend_order_tr, ~] = f_hcluster(raster', 'cosine', 1);

raster_sort_cell = raster(dend_order_cell,:);
raster_sort_ctr = raster_sort_cell(:,dend_order_tr);
figure;
imagesc(raster_sort_ctr);
title('model raster sort sort');

%% compute ensemble cells 

reliab_thresh = 0.8;
ens_size = [10 10 10];
cell_overlap_fraction = 0.2;

ens_list = cell(numel(ens_size),1);
cell_list = find(cell_reliability<=reliab_thresh);
for n_ens = 1:(numel(ens_size)-1)
    num_corr_cells = round(ens_size(n_ens)*cell_overlap_fraction);
    for n_ens2 = (n_ens+1):(numel(ens_size))
        chosen_cells = randsample(cell_list, num_corr_cells);
        ens_list{n_ens} = cat(1, ens_list{n_ens}, chosen_cells);
        ens_list{n_ens2} = cat(1, ens_list{n_ens2}, chosen_cells);
        cell_list(logical(sum(cell_list==chosen_cells',2))) = [];
    end
end
for n_ens = 1:numel(ens_size)
    chosen_cells = randsample(cell_list,ens_size(n_ens) - numel(ens_list{n_ens}));
    ens_list{n_ens} = cat(1, ens_list{n_ens}, chosen_cells);
    cell_list(logical(sum(cell_list==chosen_cells',2))) = [];
end

%% add to raster
raster_ens = raster;

num_corr_trials = [8 8 8];
available_trials = active_trials;
ens_trials = cell(numel(ens_size),1);
for n_ens = 1:numel(ens_size)
    chosen_trials = randsample(num_trials, num_corr_trials(n_ens));
    ens_trials{n_ens} = chosen_trials;
    for n_cell_ind = 1:ens_size(n_ens)
        n_cell = ens_list{n_ens}(n_cell_ind);
        for n_tr = 1:numel(chosen_trials)
            tr1 = chosen_trials(n_tr);
            if ~raster_ens(n_cell,tr1)
                if numel(available_trials{n_cell}) > 0
                    if numel(available_trials{n_cell})>1
                        chosen_tr = randsample(available_trials{n_cell}, 1);
                    else
                        chosen_tr = available_trials{n_cell};
                    end
                    raster_ens(n_cell,tr1) = raster_ens(n_cell,chosen_tr);
                    raster_ens(n_cell,chosen_tr) = 0;
                    available_trials{n_cell}(available_trials{n_cell} == chosen_tr) = [];
                else
                    raster_ens(n_cell,tr1) = randsample(peak_mag_list, 1);
                end
            else
                available_trials{n_cell}(available_trials{n_cell} == tr1) = [];
            end
        end
    end
end

%% analysis
figure;
imagesc(raster_ens)
title('model raster ens');

[dend_order_cell2, ~] = f_hcluster(raster_ens, 'cosine', 1);
[dend_order_tr2, ~] = f_hcluster(raster_ens', 'cosine', 1);

raster_ens_sort_cell = raster_ens(dend_order_cell2,:);
raster_ens_sort_ctr = raster_ens_sort_cell(:,dend_order_tr2);
figure;
imagesc(raster_ens_sort_ctr);
title('model raster ens sort sort');

ens_vals = cell(numel(ens_list),1);
for n_ens1 = 1:numel(ens_list)
    ens_vals{n_ens1} = cell(ens_size(n_ens1),1);
    for n_cell_ind1 = 1:ens_size(n_ens1)
        ens_vals{n_ens1}{n_cell_ind1} = raster_ens(ens_list{n_ens1}(n_cell_ind1),ens_trials{n_ens1});
    end
end

ens_list_sort = cell(numel(ens_size),1);
ens_trials_sort = cell(numel(ens_size),1);
for n_ens = 1:numel(ens_list)
    fprintf('ensemble %d cells', n_ens);
    ens_list_sort{n_ens} = find(sum(dend_order_cell2 == ens_list{n_ens}));
    ens_list_sort{n_ens}
    fprintf('ensemble %d trials', n_ens);
    ens_trials_sort{n_ens} = find(sum(dend_order_tr2 == ens_trials{n_ens},1));
    ens_trials_sort{n_ens}
end


%%
params.cond_name = 'shuff';
params.n_dset = 99;
params.normalize = 0;
f_ensemble_analysis_peaks2(raster_ens_sort_ctr, 170*ones*num_events, params, ops);
end

