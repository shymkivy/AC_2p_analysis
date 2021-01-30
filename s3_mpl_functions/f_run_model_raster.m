function params_out = f_run_model_raster(model_params, source_data, ens_params, script_params, ops)

plot_stuff = script_params.plot_stuff;
dim_est_only = script_params.dim_est_only;

num_cells = model_params.num_cells;
num_trials = model_params.num_trials;

%% make fake data
fake_raster_out = f_generate_fake_raster(model_params, source_data);

raster = fake_raster_out.raster;
raster_ens = fake_raster_out.raster_ens;
ens_cell_cell = fake_raster_out.ens_cell_cell;
ens_trials_cell = fake_raster_out.ens_trials_cell;
ens_pattern_rate = fake_raster_out.ens_pattern_rate;

%% quick analysis of rand

hc_params.num_clust = 1;
hc_params.plot_dist_mat = 0;
hclust_out_cell = f_hcluster_wrap(raster, hc_params);
hclust_out_tr = f_hcluster_wrap(raster', hc_params);

raster_sort_cell = raster(hclust_out_cell.dend_order,:);
raster_sort_ctr = raster_sort_cell(:,hclust_out_tr.dend_order);

if plot_stuff
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
    axis tight equal;
end

%%
ens_trials1 = zeros(num_trials,1);
for n_tr = 1:numel(ens_trials_cell)
    ens_trials1(ens_trials_cell{n_tr}) = n_tr;
end
ens_cells1 = zeros(num_cells,1);
for n_cell = 1:numel(ens_cell_cell)
    ens_cells1(ens_cell_cell{n_cell}) = n_cell;
end

if isempty(cat(1,ens_trials_cell{:}))
    ens_trials_cell_core = {find(ens_trials1 == 0)};
    ens_cells_cell_core = {find(ens_cells1 == 0)};
else
    ens_trials_cell_core = [{find(ens_trials1 == 0)}; ens_trials_cell];
    ens_cells_cell_core = [{find(ens_cells1 == 0)}; ens_cell_cell];
end
%% quick analysis 2

norm1 = 1;
if norm1
    raster_clust = raster_ens-mean(raster_ens,2);
else
    raster_clust = raster_ens;
end

hclust_out_cell2 = f_hcluster_wrap(raster, hc_params);
hclust_out_tr2 = f_hcluster_wrap(raster', hc_params);

raster_ens_sort_cell = raster_ens(hclust_out_cell2.dend_order,:);
raster_ens_sort_ctr = raster_ens_sort_cell(:,hclust_out_tr2.dend_order);

ens_pattern_rate_sort_cell = ens_pattern_rate(hclust_out_cell2.dend_order,:);
ens_pattern_rate_sort_ctr = ens_pattern_rate_sort_cell(:,hclust_out_tr2.dend_order);

% % spike mags
% ens_vals = cell(numel(ens_list_cell),1);
% for n_ens1 = 1:numel(ens_list_cell)
%     ens_vals{n_ens1} = cell(ens_size(n_ens1),1);
%     for n_cell_ind1 = 1:ens_size(n_ens1)
%         ens_vals{n_ens1}{n_cell_ind1} = raster_ens(ens_list_cell{n_ens1}(n_cell_ind1),ens_trials_cell{n_ens1});
%     end
% end

ens_cells_sort1 = ens_cells1(hclust_out_cell2.dend_order);
ens_trials_sort1 = ens_trials1(hclust_out_tr2.dend_order);

% ens_list_sort_cell = cell(numel(ens_size),1);
% ens_trials_sort_cell = cell(numel(ens_size),1);

%% plot

if plot_stuff
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
end
%% run real ens analysis 
% raster_ens_an = raster_ens_sort_ctr;
% ens_cells2 = ens_cells_sort1;
% ens_trials2 = ens_trials_sort1;

raster_ens_an = raster_ens;
ens_cells2 = ens_cells1;
ens_trials2 = ens_trials1;

% 
% for n_ens = 1:numel(ens_cell_cell)
%     fprintf('ensemble %d cells:  ', n_ens);
%     cells1 = find(ens_cells2 == n_ens);
%     for ii = 1:(numel(cells1))
%         fprintf('%d, ', cells1(ii));
%     end
%     fprintf('\n');
%     fprintf('ensemble %d trials: ', n_ens);
%     trials1 = find(ens_trials2 == n_ens);
%     for ii = 1:numel(trials1)
%         fprintf('%d ', trials1(ii));
%     end
%     fprintf('\n');
% end
table_data = cell(numel(ens_cell_cell),2);
for n_ens = 1:numel(ens_cell_cell)
    table_data{n_ens,1} = num2str(find(ens_cells2 == n_ens)');
    table_data{n_ens,2} = num2str(find(ens_trials2 == n_ens)');
end

if plot_stuff
    figure;
    uitable('Data', table_data , 'ColumnName', {'Ensemble cells', 'trials'}, 'Position', [20 20 500 350], 'ColumnWidth', {200});
end

if dim_est_only
    ens_out.data_dim_est = f_ensemble_comp_data_dim2(raster_ens_an, ens_params);
else
    ens_out = f_ensemble_analysis_peaks3(raster_ens_an, ens_params);
end
params_out.dimensionality_corr = ens_out.data_dim_est.dimensionality_corr;

%% evaluate ens
if ~dim_est_only
    params2.plot_stuff = plot_stuff;
    eval_out = f_evaluate_ens(ens_out, ens_trials_cell_core, ens_cells_cell_core, params2);
    params_out.cell_eval_accuracy = eval_out.clust_eval_cell.accuracy;
    params_out.trial_eval_accuracy = eval_out.clust_eval_tr.accuracy;

    %% plot 

    raster_ens_sort_cell2 = raster_ens_an(ens_out.cells.dend_order,:);
    raster_ens_sort_ctr2 = raster_ens_sort_cell2(:,ens_out.trials.dend_order);

    clust_ident_cell_align = eval_out.clust_eval_cell.aligned_seq(ens_out.cells.clust_ident+1,1);
    clust_ident_trial_align = eval_out.clust_eval_tr.aligned_seq(ens_out.trials.clust_ident+1,1);

    if plot_stuff
        figure(f1);
        subplot(4,3,10);
        imagesc(raster_ens_sort_ctr2)
        title('Post ens analysis ground truth label');
        f_plot_cell_indicator(raster_ens, ens_cells2(ens_out.cells.dend_order), ops);
        f_plot_trial_indicator2(raster, ens_trials2(ens_out.trials.dend_order), 1, ops);

        figure(f1);
        subplot(4,3,11);
        imagesc(raster_ens_sort_ctr2);
        title('Post ens analysis discovered label');
        f_plot_cell_indicator(raster_ens, clust_ident_cell_align(ens_out.cells.dend_order), ops);
        f_plot_trial_indicator2(raster, clust_ident_trial_align(ens_out.trials.dend_order), 1, ops);

        figure(f1);
        subplot(4,3,12);
        imagesc(corr(raster_ens_sort_ctr2'));
        title('Post ens analysis cell corr');
        axis tight equal
        f_plot_cell_indicator(corr(raster_ens_sort_ctr'), clust_ident_cell_align(ens_out.cells.dend_order), ops);
    end
end
end