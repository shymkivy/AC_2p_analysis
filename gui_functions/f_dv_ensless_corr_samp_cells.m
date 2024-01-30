function f_dv_ensless_corr_samp_cells(app)

disp('Computing ensless corr with cell sampling...');

sort_type = app.SorttypeDropDown.Value; %'trials', 'cells'
corr_type = app.CorrtypeDropDown.Value; % 'SI_cosine', 'SI_correlation', 'pca_dim'

equalize_tn_cell_samp = app.EqualizetrialsampsizeCheckBox.Value;

hc_params.plot_dist_mat = 0;
hc_params.plot_clusters = 0;
hc_params.num_clust = 1;

resp_type = app.ResposivecellstypeDropDown.Value;
resp_selection = app.ResponsivecellsselectDropDown.Value;

samp_range_min = app.samprangeminEditField.Value;
samp_range_max = app.samprangemaxEditField.Value;
samp_interv = app.samprangeintervEditField.Value;
num_samp = app.sampnumEditField.Value;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = size(data,1);
params = f_dv_gather_params(app);

tn_all = f_dv_get_trial_number(app);
[num_gr, num_tn] = size(tn_all);
trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);

if equalize_tn_cell_samp % sampling from cells from every tn separately, in equal amounts, and then combining.. limited by red tn size
    samp_range_min = round(samp_range_min/num_tn);
    samp_range_max = round(samp_range_max/num_tn);
    samp_interv = round(samp_interv/num_tn);
end

corr_vals = cell(num_dsets, num_gr, num_tn);
num_resp_cells = zeros(num_dsets, num_gr, num_tn);
max_range = 0;
fprintf('Dsets n/%d: ', num_dsets);
for n_dset = 1:num_dsets
    fprintf('%d..', n_dset);
    ddata = data(n_dset,:);
    [cdata, stats1] = f_dv_get_new_cdata_stats(app, ddata, params);
    params.cdata = cdata;
    mmn_freq = ddata.MMN_freq{1};

    firing_rate = cat(1,cdata.S_sm);
    
    for n_gr = 1:num_gr
        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all(n_gr,:));
       
        resp_cells_split = f_dv_get_resp_vals_cells(app, stats1, tn_all(n_gr,:), resp_type, 'resp split');
        
        trial_types = ddata.trial_types{1};
        stim_frame_index = ddata.stim_frame_index{1};
        [~, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);
        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_frames);
        
        if ~isempty(mmn_freq)
            trial_types_ctx2 = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);
            trial_types_all = [trial_types, trial_types_ctx2];
        else
            trial_types_all = trial_types;
        end
    
        if equalize_tn_cell_samp
            num_cells_tn = sum(resp_cells_split,1);
            min_num_cells = min(num_cells_tn);
            [resp_cell_split2,tn_idx] = find(resp_cells_split);
        end
    
        for n_tn = 1:num_tn
    
            tr_idx = logical(sum(app.ops.context_types_all(tn_all(n_gr,n_tn)) == trial_types_all,2));
            tr_data = trial_data_sort(:,:,tr_idx);
            
            resp_cells2 = find(resp_cells(:, n_tn));
            
            if ~equalize_tn_cell_samp
                min_num_cells = sum(resp_cells(:, n_tn));
            end
    
            [num_cells2, ~, num_tr] = size(tr_data);
              
            %firing_rate2 = firing_rate(resp_cells1,:);
            if min_num_cells >=2
                
                temp_max = min([samp_range_max, min_num_cells]);
                samp_range = samp_range_min:samp_interv:temp_max;
                max_range = max([temp_max, max_range]);
                
                num_range = numel(samp_range);
                
                SI_all = zeros(num_range, num_samp);
                for n_range = 1:num_range
                    for n_samp = 1:num_samp
                        % add step where we equalize the number of cells taken
                        % from each
                        
                        if equalize_tn_cell_samp
                            samp_cell2 = cell(num_tn,1);
                            for n_tn2 = 1:num_tn
                                resp_cell_split3 = find(resp_cells_split(:,n_tn2));
                                samp_cell2{n_tn2} = randsample(resp_cell_split3, samp_range(n_range));
                            end
                            samp_cells = unique(cat(1, samp_cell2{:}));
                        else
                            samp_cells = randsample(resp_cells2, samp_range(n_range));
                        end
                        
                        tr_data2 = tr_data(samp_cells,:,:);
                        
                        if strcmpi(sort_type, 'cells')
                            tr_data2_2d = reshape(tr_data2, samp_range(n_range), []);
                        elseif strcmpi(sort_type, 'trials')
                            tr_data2_2d = reshape(tr_data2, [], num_tr)';
                        end
                        
                        if strcmpi(corr_type, 'SI_cosine')
                            SI_vals = 1 - pdist(tr_data2_2d, 'cosine');
                            SI_out = mean(SI_vals);
                        elseif strcmpi(corr_type, 'SI_correlation')
                            SI_vals = 1 - pdist(tr_data2_2d, 'correlation');
                            SI_out = mean(SI_vals);
                        elseif strcmpi(corr_type, 'pca_dim')
                            dim_est_pca = f_ensemble_comp_data_dim2(tr_data2_2d, params.est_params_pca);
                            SI_out = dim_est_pca.dimensionality_corr;
                        end
    
                        SI_all(n_range, n_samp) = SI_out;
                    end
                end
                
                num_resp_cells(n_dset, n_gr, n_tn) = num_cells2;
                corr_vals{n_dset, n_gr, n_tn} = mean(SI_all,2);
            end
        end
    end
end
fprintf('\n');

if equalize_tn_cell_samp
    temp_x = (samp_range_min:samp_interv:max_range)*num_tn;
else
    temp_x = samp_range_min:samp_interv:max_range;
end


num_x = numel(temp_x);
dset_mean = nan(num_tn, num_x);
dset_sem = nan(num_tn, num_x);
for n_tn = 1:num_tn
    corr_vals2 = nan(num_dsets, num_gr, num_x);
    for n_gr = 1:num_gr
        for n_dset = 1:num_dsets
            num_val = numel(corr_vals{n_dset, n_gr, n_tn});
            corr_vals2(n_dset, n_gr, 1:num_val) = corr_vals{n_dset, n_gr, n_tn};
        end
    end
    corr_vals3 = reshape(corr_vals2, [num_dsets*num_gr, num_x]);
    for n_x = 1:num_x
        idx1 = ~isnan(corr_vals3(:, n_x));
        if sum(idx1)
            temp_vals = corr_vals3(idx1, n_x);
            dset_mean(n_tn, n_x) = mean(temp_vals);
            dset_sem(n_tn, n_x) = std(temp_vals)/sqrt(max(numel(temp_vals) - 1, 1));
        end
    end
end

x_lab_categories = categorical(app.ops.context_types_labels(tn_all));
color1 = app.ops.context_types_all_colors2;

figure; hold on;
for n_tn = 1:num_tn
    temp_mean = dset_mean(n_tn, :);
    temp_sem = dset_sem(n_tn, :);
    plot(temp_x, temp_mean, '-', 'color', color1{tn_all(1, n_tn)}, 'linewidth', 2);
    errorbar(temp_x, temp_mean, temp_sem, 'color', color1{tn_all(1, n_tn)})
end
axis tight;
ylabel(corr_type, 'interpreter', 'none'); xlabel('cells');
title(sprintf('Corr samp; sort by %s, corr %s; %s; %s; %s', sort_type, corr_type, title_tag, resp_type, resp_selection), 'interpreter', 'none');

figure; hold on;
for n_tn = 1:num_tn
    temp_mean = dset_mean(n_tn, :);
    temp_sem = dset_sem(n_tn, :);
    shadedErrorBar_YS(temp_x, temp_mean, temp_sem, color1{tn_all(1, n_tn)});
end
axis tight;
ylabel(corr_type, 'interpreter', 'none'); xlabel('cells');
title(sprintf('Corr samp; sort by %s, corr %s; %s; %s; %s', sort_type, corr_type, title_tag, resp_type, resp_selection), 'interpreter', 'none');

% 
% figure; hold on;
% for n_tn = 1:numel(tn_all)
%     temp_x = samp_range_min:samp_interv:max_range;
%     temp_mean = mean(corr_vals{1, n_tn}, 2);
%     temp_sem = std(corr_vals{1, n_tn}, [], 2)/sqrt(num_samp-1);
%     plot(temp_x, temp_mean, '-', 'color', color1{tn_all(n_tn)}, 'linewidth', 2)
%     errorbar(temp_x, temp_mean, temp_sem, 'color', color1{tn_all(n_tn)})
% end
% 
% 
% figure; hold on;
% for n_tn = 1:numel(tn_all)
%     plot(num_resp_cells(:,n_tn), corr_vals(:,n_tn), 'o', 'color', color1{tn_all(n_tn)}, 'linewidth', 2)
% end
% xlabel('num cells'); ylabel('Pairwise correlation')
% title(sprintf('Mean pairwise correlations; %s', title_tag), 'interpreter', 'none');
% legend(x_lab_categories)
% 
% vals1 = corr_vals;
% figure; hold on;
% bar(categorical(x_lab_categories,x_lab_categories), mean(vals1, 1));
% for n_tn = 1:numel(tn_all)
%     bar(x_lab_categories(n_tn), mean(vals1(:,n_tn), 1), 'FaceColor', min(color1{tn_all(n_tn)}+.5, 1), 'EdgeColor', color1{tn_all(n_tn)}, 'linewidth', 2)
%     errorbar(x_lab_categories(n_tn), mean(vals1(:,n_tn), 1), std(vals1(:,n_tn), [], 1)/sqrt(num_dsets-1), 'k.', 'linewidth', 2)
% end
% ylabel('pairwise correlation');
% title(sprintf('Mean correlation across datasets, %s', title_tag), 'interpreter', 'none')
% 
% vals1 = corr_vals./sum(num_resp_cells, 1);
% figure; hold on;
% bar(categorical(x_lab_categories,x_lab_categories), mean(vals1, 1));
% for n_tn = 1:numel(tn_all)
%     bar(x_lab_categories(n_tn), mean(vals1(:,n_tn), 1), 'FaceColor', min(color1{tn_all(n_tn)}+.5, 1), 'EdgeColor', color1{tn_all(n_tn)}, 'linewidth', 2)
%     errorbar(x_lab_categories(n_tn), mean(vals1(:,n_tn), 1), std(vals1(:,n_tn), [], 1)/sqrt(num_dsets-1), 'k.', 'linewidth', 2)
% end
% ylabel('pairwise correlation / num resp cells');
% title(sprintf('Mean correlation per cells across datasets, %s', title_tag), 'interpreter', 'none')
% 
% vals2 = sum(num_resp_cells, 1);
% figure; hold on;
% bar(categorical(x_lab_categories,x_lab_categories), vals2);
% for n_tn = 1:numel(tn_all)
%     bar(x_lab_categories(n_tn), vals2(n_tn), 'FaceColor', min(color1{tn_all(n_tn)}+.5, 1), 'EdgeColor', color1{tn_all(n_tn)}, 'linewidth', 2)
% end
% ylabel('num responsive cells');
% title(sprintf('Numebr of cells across datasets, %s', title_tag), 'interpreter', 'none')
% 
% 
% % figure; hold on
% % for n_freq = 1:numel(tn_all)
% %     plot([0.5 1 2 4], corr_vals(:,n_freq), 'o-', 'color', color1{tn_all(n_freq), 'linewidth', 2)
% % end
% % xlabel('ISI duration'); ylabel('Pairwise correlation')
% % title('Mean pairwise correlations');
% % 
% % figure; imagesc(reshape(color1, [10 1 3]))

disp('Done');

% sp = 1:4:21;
% x = gca;
% for n_sp = 1:numel(sp)
%     x.Children(sp(n_sp)).LineWidth = 1;
% end


end