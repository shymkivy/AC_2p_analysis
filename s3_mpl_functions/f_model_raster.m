function f_model_raster(data, ops)

num_reps = 10;

%% ensemble params
model_params.num_cells = [40 50 60];
model_params.num_trials = 50;

model_params.reliab_thresh = 0.8;
model_params.ens_size = {[5 5] [5 5 5]};
model_params.cell_overlap_fraction = 0;
model_params.num_corr_trials = {[6 6] [6 6 6 ]};

plot_stuff = 0;

%% clust params
ens_params.cond_name = 'fake data';
ens_params.n_dset = 99;
ens_params.normalize = 'norm_full'; %'norm_full' 'norm_mean' 'none'
ens_params.num_comps = [];
ens_params.plot_stuff = 0;
ens_params.ensamble_method = 'nmf';
ens_params.use_LR_proj = 0;
ens_params.ensamble_extraction = 'clust'; % 'thresh', 'clust'
ens_params.ensamble_extraction_thresh = 'signal_clust_thresh'; % 'shuff'. 'clust_thresh', 'signal_z'

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

source_data.peak_mag_list = cat(1,peak_mag_list{:});
% figure; histogram(peak_mag_list, 20);
% [f, x] = ksdensity(peak_mag_list);
% figure;
% plot(x,f)
source_data.reliab_list = cat(1,reliab_list{:});
% x = sort(reliab_list);
% f = (1:numel(reliab_list))/numel(reliab_list);
% figure; plot(x, f)

%% build input struct
input_params = struct();
pindx = 1;
for n_cell = 1:numel(model_params.num_cells)
    for n_tr = 1:numel(model_params.num_trials)
        for n_ens = 1:numel(model_params.ens_size)
            for n_rep = 1:num_reps
                input_params(pindx).num_cells = model_params.num_cells(n_cell);
                input_params(pindx).num_trials = model_params.num_trials(n_tr);
                input_params(pindx).n_rep = n_rep;
                input_params(pindx).reliab_thresh = model_params.reliab_thresh;
                input_params(pindx).ens_size = model_params.ens_size{n_ens};
                input_params(pindx).cell_overlap_fraction = model_params.cell_overlap_fraction;
                input_params(pindx).num_corr_trials = model_params.num_corr_trials{n_ens};
                input_params(pindx).num_ens = numel(model_params.ens_size{n_ens});
                pindx = pindx +1;
            end
        end
    end
end
%% run model
eval_out = input_params;
for n_rep = 1:numel(input_params)
    x = f_run_model_raster(input_params(n_rep), source_data, ens_params, ops, plot_stuff);
    eval_out(n_rep).eval_num_clust = x.eval_num_clust;
    eval_out(n_rep).cell_eval_accuracy = x.cell_eval_accuracy;
    eval_out(n_rep).trial_eval_accuracy = x.trial_eval_accuracy;
end
%%


%%
disp('Done')
end

