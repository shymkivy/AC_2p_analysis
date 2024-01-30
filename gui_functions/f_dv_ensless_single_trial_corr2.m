function f_dv_ensless_single_trial_corr2(app)

disp('Computing single trial corr...');

hc_params.plot_dist_mat = 0;
hc_params.plot_clusters = 0;
hc_params.num_clust = 1;

dist_method = lower(app.DistmethodDropDown.Value);
do_sim = app.dosimilarityCheckBox.Value;
dist_ref = app.DistreferenceDropDown.Value;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = size(data,1);
params = f_dv_gather_params(app);

tn_all = f_dv_get_trial_number(app);
[num_gr, num_tn] = size(tn_all);
trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);

corr_vals = zeros(num_dsets, num_gr, num_tn);
num_resp_cells = zeros(num_dsets, num_gr, num_tn);
has_data = false(num_dsets, num_gr, num_tn);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    [cdata, stats1] = f_dv_get_new_cdata_stats(app, ddata, params);
    mmn_freq = ddata.MMN_freq{1};
    firing_rate = cat(1,cdata.S_sm);
  
    
    
    trial_types = ddata.trial_types{1};
    stim_frame_index = ddata.stim_frame_index{1};
    [plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);
    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_frames);
    
    if ~isempty(mmn_freq)
        trial_types_ctx2 = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);
        trial_types_all = [trial_types, trial_types_ctx2];
    else
        trial_types_all = trial_types;
    end
    
    for n_gr = 1:num_gr
        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all(n_gr,:));
        resp_cells_split = f_dv_get_resp_vals_cells(app, stats1, tn_all(n_gr,:), [], 'resp split');
        for n_tn = 1:num_tn
            tr_idx = logical(sum(app.ops.context_types_all(tn_all(n_gr, n_tn)) == trial_types_all,2));
            resp_cells1 = resp_cells(:, n_tn);
    
            tr_data = trial_data_sort(:,:,tr_idx);
            tr_data2 = tr_data(resp_cells1,:,:);
            
            %tr_data_2d_tr = reshape(tr_data2, [], num_tr);
            tr_data_2d_tr = squeeze(mean(tr_data2,2));
    
            firing_rate2 = firing_rate(resp_cells1,:);
            [num_cells2, ~, num_tr] = size(tr_data2);
            
            if num_cells2>1
                if strcmpi(dist_ref, 'pairwise')
                    corr_vals1 = pdist(tr_data_2d_tr', dist_method);
                else
                    if strcmpi(dist_ref, 'zero')
                        ref_vec = zeros(num_cells2, 1);
                    elseif strcmpi(dist_ref, 'trial ave')
                        ref_vec = mean(tr_data_2d_tr,2);
                    end
                    corr_vals1 = pdist2(ref_vec', tr_data_2d_tr', dist_method);
                end
                num_resp_cells(n_dset, n_gr, n_tn) = num_cells2;
                num_resp_cells_split(n_dset, n_gr, n_tn) = sum(resp_cells_split(:,n_tn));
                corr_vals(n_dset, n_gr, n_tn) = mean(corr_vals1);
                has_data(n_dset, n_gr, n_tn) = 1;
            end
        end
    end
end

corr_vals2 = corr_vals;
ytag1 = '';
if strcmpi(dist_method, 'euclidean')
    if app.normalizeeucCheckBox.Value
        corr_vals2(has_data) = corr_vals(has_data)./sqrt(num_resp_cells(has_data))*mean(sqrt(num_resp_cells(has_data)));
        ytag1 = 'normalized ';
    end
end

if do_sim
    corr_vals2 = 1 - corr_vals2;
    ytag2 = ' similarity';
else
    ytag2 = ' distance';
end

if ~strcmpi(dist_ref, 'pairwise')
    ytag4 = '';
    ytag3 = sprintf(' from %s', dist_ref);
else
    ytag4 = 'pairwise ';
    ytag3 = '';
end

ytag = sprintf('%s%s%s%s%s', ytag4, ytag1, dist_method, ytag2, ytag3);

title_tag2 = sprintf('%s; %s; %s', title_tag, app.ResposivecellstypeDropDown.Value, app.ResponsivecellsselectDropDown.Value);

x_lab_categories = categorical(app.ops.context_types_labels(tn_all(1,:)));
color1 = app.ops.context_types_all_colors2;
pl = cell(num_tn,1);
figure; hold on;
for n_tn = 1:num_tn
    has_data2 = has_data(:,:,n_tn);
    num_cells2 = num_resp_cells_split(:,:,n_tn);
    corr_vals3 = corr_vals2(:, :, n_tn);
    pl{n_tn} = plot(num_cells2(has_data2), corr_vals3(has_data2), '.', 'color', color1{tn_all(1,n_tn)}, markerSize=20);
    plot(num_cells2(has_data2), corr_vals3(has_data2), 'o', color='k', markerSize=5, lineWidth=0.5);
end
xlabel('num context cells'); ylabel(ytag)
title(sprintf('%s\n%s vs cell count', title_tag2, ytag), 'interpreter', 'none');
legend([pl{:}], x_lab_categories)

y_data = cell(num_tn, 1);
for n_tn = 1:num_tn
    has_data2 = has_data(:,:,n_tn);
    corr_vals3 = corr_vals2(:, :, n_tn);
    y_data{n_tn} = corr_vals3(has_data2);
end
cont_idx = logical(sum(tn_all==[18;28],1));
dd_idx = logical(sum(tn_all==[20;30],1));
if sum(dd_idx) && sum(cont_idx) && num_dsets>1
    p_rs= ranksum(y_data{cont_idx},y_data{dd_idx});
    [~, p_t] = ttest2(y_data{cont_idx},y_data{dd_idx});
    stats_tag = sprintf('; prs=%.4f; pt=%.4f', p_rs, p_t);
else
    stats_tag = '';
end

figure; hold on;
bar(categorical(x_lab_categories,x_lab_categories), zeros(num_tn, 1));
for n_tn = 1:num_tn
    mean1 = mean(y_data{n_tn});
    sem1 = std(y_data{n_tn})/sqrt(max(num_dsets-1,1));
    bar(x_lab_categories(n_tn), mean1, 'FaceColor', min(color1{tn_all(1,n_tn)}+.5, 1), 'EdgeColor', color1{tn_all(1,n_tn)}, 'linewidth', 2)
    errorbar(x_lab_categories(n_tn), mean1, sem1, 'k.', 'linewidth', 2)
end
ylabel(ytag);
title(sprintf('%s\n%s%s', title_tag2, ytag, stats_tag), 'interpreter', 'none')

figure; hold on;
bar(categorical(x_lab_categories,x_lab_categories), zeros(num_tn, 1));
for n_tn = 1:num_tn
    mean1 = mean(y_data{n_tn});
    std1 = std(y_data{n_tn});
    sem1 = std1/sqrt(max(num_dsets-1,1));
    num_dsets2 = numel(y_data{n_tn});
    x_data = (rand(num_dsets2,1)-0.5)/5+n_tn;
    scatter(x_data, y_data{n_tn}, '.', MarkerEdgecolor=color1{tn_all(1,n_tn)});
    %plot(x_data, corr_vals3(has_data2), 'o', color='k', markerSize=5, lineWidth=0.5)
    %bar(x_lab_categories(n_tn), mean1, EdgeColor='k', linewidth=2, faceColor='none'); %color1{tn_all(1,n_tn)}
    plot(n_tn, mean1, '_', color='k', linewidth=2, markersize=15);
    errorbar(x_lab_categories(n_tn), mean1, sem1, 'k.', linewidth=1, markersize=10);
end
ylabel(ytag);
title(sprintf('%s\n%s 2%s', title_tag2, ytag, stats_tag), 'interpreter', 'none')

figure; hold on;
bar(categorical(x_lab_categories,x_lab_categories), zeros(num_tn, 1));
for n_tn = 1:num_tn
    num_cells2 = num_resp_cells_split(:,:,n_tn);
    bar(x_lab_categories(n_tn), sum(num_cells2(:)), 'FaceColor', min(color1{tn_all(1,n_tn)}+.5, 1), 'EdgeColor', color1{tn_all(1,n_tn)}, 'linewidth', 2)
end
ylabel('num responsive cells');
title(sprintf('%s\n%s', title_tag2, ytag), 'interpreter', 'none')


% figure; hold on
% for n_freq = 1:numel(tn_all)
%     plot([0.5 1 2 4], corr_vals(:,n_freq), 'o-', 'color', color1{tn_all(n_freq), 'linewidth', 2)
% end
% xlabel('ISI duration'); ylabel('Pairwise correlation')
% title('Mean pairwise correlations');
% 
% figure; imagesc(reshape(color1, [10 1 3]))

disp('Done')

% testing how corr increasese with pop size
% 
% pop_all = 1:5:500;
% 
% dist_all = zeros(numel(pop_all),1);
% 
% for n_pop = 1:numel(pop_all)
% 
%     vecd = rand(pop_all(n_pop),1);
% 
%     vec0 = zeros(pop_all(n_pop),1);
%     vec1 = ones(pop_all(n_pop),1);
% 
%     %dist_all(n_pop) = sqrt(mean((vecd).^2));
% 
%     %dist_all(n_pop) = norm(vecd, 'fro')/norm(vec1, 'fro');
% 
%     dist_all(n_pop) = pdist2(vec1', vecd', 'euclidean')/pdist2(vec1', vec0', 'euclidean');
% end
% 
% 
% figure;
% plot(pop_all, dist_all)

end