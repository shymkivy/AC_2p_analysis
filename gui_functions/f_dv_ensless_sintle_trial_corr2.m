function f_dv_ensless_sintle_trial_corr2(app)

disp('Computing single trial corr...');

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = size(data,1);

params = f_dv_gather_params(app);

tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';

select_resp_cells = app.selectrespcellsCheckBox.Value;
resort_by_ens = app.resortbyensCheckBox.Value;
sort_trials = app.sorttrialsCheckBox.Value;
sort_with_full_firing_rate = app.sortwithfullfrCheckBox.Value;

corr_vals = zeros(num_dsets, numel(tn_all));
num_cells_all = zeros(num_dsets, numel(tn_all));
num_resp_cells = zeros(num_dsets, numel(tn_all));
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    
    if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
        cdata = f_dv_compute_cdata(ddata, params);
        stats1 = ddata.stats{app.mplSpinner.Value};
    else
        cdata = cell(ddata.num_planes,1);
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            cdata{n_pl} = f_dv_compute_cdata(ddata, params);
        end
        cdata = cat(1,cdata{:});
        stats1 = cat(1,ddata.stats{:});
    end
    
    firing_rate = cat(1,cdata.S_sm);
    num_cells = sum([cdata.num_cells]);
    
    [~, resp_cells_all] = f_dv_get_resp_vals_cells(stats1, tn_all, app.ResposivecellstypeDropDown.Value, app.LimitresptrialsCheckBox.Value, app.RespthreshEditField.Value);
   
    for n_tn = 1:numel(tn_all)
        
        
        
%         if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
%             resp_cell = logical(sum(stats1{n_pl}.resp_cells_peak(:,tn_all),2));
%         else
%             resp_cells_all = cell(numel(stats1),1);
%             for n_pl2 = 1:numel(stats1)
%                 resp_cells_all{n_pl2} = logical(sum(stats1{n_pl2}.resp_cells_peak(:,tn_all),2));
%             end
%             resp_cell = cat(1,resp_cells_all{:});
%         end
        
        %%
        trial_types = ddata.trial_types{1};
        stim_frame_index = ddata.stim_frame_index{1};
        
        trial_window = f_str_to_array(app.plot_BaserespwinEditField.Value);
        [plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_frames);
        if ~isempty(ddata.MMN_freq{1})
            [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, ddata.MMN_freq{1}, app.ops);
        else
            trial_data_sort_wctx = trial_data_sort;
            trial_types_wctx = trial_types;
        end

        tr_idx = trial_types_wctx == app.ops.context_types_all(tn_all(n_tn));
        tr_data = trial_data_sort_wctx(:,:,tr_idx);
        
        num_resp_cells(n_dset, n_tn) = sum(resp_cells_all{n_tn});
        if select_resp_cells
            resp_cell = resp_cells_all{n_tn};
            tr_data2 = tr_data(resp_cell,:,:);
            firing_rate2 = firing_rate(resp_cell,:);
        else
            tr_data2 = tr_data;
            firing_rate2 = firing_rate;
        end

        [num_cells2, ~, num_tr] = size(tr_data2);
        num_cells_all(n_dset, n_tn) = num_cells2;
        
        hc_params.plot_dist_mat = 0;
        hc_params.plot_clusters = 0;
        hc_params.num_clust = 1;

        tr_data_2d_tr = reshape(tr_data2, [], num_tr);
        hclust_out_trial = f_hcluster_wrap(tr_data_2d_tr', hc_params);
        tr_data3 = tr_data2(:,:,hclust_out_trial.dend_order);

        SI = 1-hclust_out_trial.dist;
        SI_vals = tril(SI,-1);
        SI_vals(SI_vals==0) = [];
        
        corr_vals(n_dset, n_tn) = mean(SI_vals);
   
    end
end

x_lab_categories = categorical(app.ops.context_types_labels(tn_all));
color1 = app.ops.context_types_all_colors2;
figure; hold on;
for n_freq = 1:numel(tn_all)
    plot(num_cells_all(:,n_freq), corr_vals(:,n_freq), 'o', 'color', color1{tn_all(n_freq)}, 'linewidth', 2)
end
xlabel('num cells'); ylabel('Pairwise correlation')
title(sprintf('Mean pairwise correlations; %s', title_tag), 'interpreter', 'none');
legend(x_lab_categories)

vals1 = corr_vals;
figure; hold on;
bar(categorical(x_lab_categories,x_lab_categories), mean(vals1));
for n_freq = 1:numel(tn_all)
    bar(x_lab_categories(n_freq), mean(vals1(:,n_freq)), 'FaceColor', min(color1{tn_all(n_freq)}+.5, 1), 'EdgeColor', color1{tn_all(n_freq)}, 'linewidth', 2)
    errorbar(x_lab_categories(n_freq), mean(vals1(:,n_freq)), std(vals1(:,n_freq))/sqrt(num_dsets-1), 'k.', 'linewidth', 2)
end
ylabel('pairwise correlation');
title(sprintf('Mean correlation across datasets, %s', title_tag), 'interpreter', 'none')

vals1 = corr_vals./sum(num_resp_cells);
figure; hold on;
bar(categorical(x_lab_categories,x_lab_categories), mean(vals1));
for n_freq = 1:numel(tn_all)
    bar(x_lab_categories(n_freq), mean(vals1(:,n_freq)), 'FaceColor', min(color1{tn_all(n_freq)}+.5, 1), 'EdgeColor', color1{tn_all(n_freq)}, 'linewidth', 2)
    errorbar(x_lab_categories(n_freq), mean(vals1(:,n_freq)), std(vals1(:,n_freq))/sqrt(num_dsets-1), 'k.', 'linewidth', 2)
end
ylabel('pairwise correlation / num resp cells');
title(sprintf('Mean correlation per cells across datasets, %s', title_tag), 'interpreter', 'none')

vals2 = sum(num_resp_cells);
figure; hold on;
bar(categorical(x_lab_categories,x_lab_categories), vals2);
for n_freq = 1:numel(tn_all)
    bar(x_lab_categories(n_freq), vals2(n_freq), 'FaceColor', min(color1{tn_all(n_freq)}+.5, 1), 'EdgeColor', color1{tn_all(n_freq)}, 'linewidth', 2)
end
ylabel('num responsive cells');
title(sprintf('Numebr of cells across datasets, %s', title_tag), 'interpreter', 'none')


% figure; hold on
% for n_freq = 1:numel(tn_all)
%     plot([0.5 1 2 4], corr_vals(:,n_freq), 'o-', 'color', color1{tn_all(n_freq), 'linewidth', 2)
% end
% xlabel('ISI duration'); ylabel('Pairwise correlation')
% title('Mean pairwise correlations');
% 
% figure; imagesc(reshape(color1, [10 1 3]))

disp('Done')

end