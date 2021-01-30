function f_model_raster(data, ops)

script_params.num_reps = 1;
script_params.plot_stuff = 1;
script_params.dim_est_only = 0;

%% ensemble params
model_params.num_cells = 40;
model_params.num_trials = 50;

model_params.reliab_thresh = 0.8;
model_params.ens_size = {[5 5 5]};
model_params.cell_overlap_fraction = 0;
model_params.num_corr_trials = {[6 6 6]};

%% clust params
ens_params.cond_name = 'fake data';
ens_params.n_dset = 99;
ens_params.normalize = 'norm_full'; %'norm_full' 'norm_mean' 'none'
ens_params.corr_comp_thresh = .85;
ens_params.num_comps = [];
ens_params.plot_stuff = 0;
ens_params.ensamble_method = 'nmf';
ens_params.use_LR_proj = 0;
ens_params.ensamble_extraction = 'thresh'; % 'thresh', 'clust'
ens_params.ensamble_extraction_thresh = 'signal_clust_thresh'; % 'shuff'. 'clust_thresh', 'signal_z'

%load('reliability_dd12')

%% load cell distributions

tt = [20, 30];

reliab_list = cell(numel(ops.regions_to_analyze), 1);
peak_mag_list = cell(numel(ops.regions_to_analyze), 1);
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    
    reliab_list1 = cell(numel(cdata.area), numel(tt));
    peak_mag_list1 = cell(numel(cdata.area), numel(tt));
    for n_dset = 1:numel(cdata.area)
        trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
        trial_types = cdata.trial_types_wctx{n_dset};
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
            trial_peaks3 = trial_peaks2(trial_peaks2>cell_reliability_thr2);
            peak_mag_list1{n_dset, n_tt} = trial_peaks3(:);
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
            for n_rep = 1:script_params.num_reps
                input_params(pindx).num_cells = model_params.num_cells(n_cell);
                input_params(pindx).num_trials = model_params.num_trials(n_tr);
                input_params(pindx).n_rep = n_rep;
                input_params(pindx).reliab_thresh = model_params.reliab_thresh;
                input_params(pindx).ens_size = model_params.ens_size{n_ens};
                input_params(pindx).cell_overlap_fraction = model_params.cell_overlap_fraction;
                input_params(pindx).num_corr_trials = model_params.num_corr_trials{n_ens};
                input_params(pindx).num_ens = sum(model_params.ens_size{n_ens}>0);
                pindx = pindx +1;
            end
        end
    end
end

%% run model
eval_out = input_params;
for n_rep = 1:numel(input_params)
    x = f_run_model_raster(input_params(n_rep), source_data, ens_params, script_params, ops);
    eval_out(n_rep).dimensionality_corr = x.dimensionality_corr;
    if ~script_params.dim_est_only
        eval_out(n_rep).cell_eval_accuracy = x.cell_eval_accuracy;
        eval_out(n_rep).trial_eval_accuracy = x.trial_eval_accuracy;
    end
end
%%

figure; hold on;
bar(mean([eval_out.dimensionality_corr]))
errorbar(1, mean([eval_out.dimensionality_corr]), std([eval_out.dimensionality_corr])/sqrt(numel([eval_out.dimensionality_corr])-1))
title(sprintf('mean corr dim = %.2f, real ens num = %.2f', mean([eval_out.dimensionality_corr]), mean([eval_out.num_ens])))
%%
disp('Done')
end

