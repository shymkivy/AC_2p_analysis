function f_dv_trial_ave_pca(app)

normalize1 = 0;
plot_extra = 0;

n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = numel(data.experiment);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

num_t = sum(trial_frames);

tn_all = f_dv_get_trial_number(app);
num_tn = numel(tn_all);

params = f_dv_gather_params(app);

resp_all = cell(num_dsets, 1);
cell_counts = zeros(num_dsets, 1);

reg_all = app.ops.regions_to_analyze;
[region_num_all, reg_tag] = f_dv_get_region_sel_val2(app);

num_groups = size(region_num_all,1);

group_data = cell(num_groups,1);
num_cells_all = zeros(num_groups,1);
for n_gr = 1:num_groups
    fprintf('group %d/%d; dset #/%d: ', n_gr, num_groups, num_dsets);
    region_num = region_num_all(n_gr,:);
    for n_dset = 1:num_dsets
        fprintf('..%d', n_dset);
        data1 =  data(n_dset,:);
        stats1 = data1.stats{n_pl};
        params.n_dset = find(data1.idx == app.data.idx);

        cdata = f_dv_compute_cdata(data1, params);

        num_cells = sum([stats1.num_cells]);

        firing_rate = cdata.S_sm;
        trial_types = data1.trial_types{1};
        stim_times = data1.stim_frame_index{n_pl};
        mmn_freq = data1.MMN_freq{1};

        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all);
        %cell_is_resp = stats1.peak_resp_cells(:,tn_all);

        if app.UseregdatalabelsCheckBox.Value
            if ~isempty(data1.registered_data{1})
                reg_cell_labels = data1.registered_data{1}.reg_labels;
            else
                reg_cell_labels = zeros(num_cells,1);
            end
        else
            reg_idx = find(strcmpi(reg_all, data1.area));
            reg_cell_labels = ones(num_cells,1)*reg_idx;
        end

        reg_cell_idx = logical(sum(reg_cell_labels == region_num,2));

        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
        [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

        % get resp cells
        resp_cell_idx = logical(sum(resp_cells,2).*reg_cell_idx);
        cell_counts(n_dset, 1) = sum(resp_cell_idx);

        resp_all{n_dset, 1} = zeros(cell_counts(n_dset, 1), num_t, num_tn);

        for n_tn = 1:num_tn
            temp_resp = trial_data_sort_wctx(resp_cell_idx,:,(trial_types_wctx == app.ops.context_types_all(tn_all(n_tn))));
            resp_all{n_dset}(:,:,n_tn) = mean(temp_resp,3);
        end
    end
    fprintf('\n');
    resp_all_pool = cat(1, resp_all{:});
    num_cells = size(resp_all_pool,1);
    
    group_data{n_gr} = resp_all_pool;
    num_cells_all(n_gr) = num_cells;
end
    
% now the pca

num_comp = 6;

dist_list = [18, 20;
             28, 30;
             2, 3;
             3, 4;
             4, 5;
             5, 6;
             6, 7;
             7, 8;
             8, 9];
         
dist_lab = {'DD-cont', 'DD-cont flip', '2-3', '3-4', '4-5', '5-6', '6-7', '7-8', '8-9'};

num_dist = size(dist_list,1);

tomp_comp_all = cell(num_groups,1);
exp_var_all = cell(num_groups,1);
for n_gr = 1:num_groups

    resp_all2d = reshape(group_data{n_gr}, num_cells_all(n_gr), []);

    if normalize1
        resp_all2d = f_normalize(resp_all2d, 'norm_mean_std');
    end
    
    %[U,S,V] = svd(resp_all2dn);
    % score*coeff'
    [~,score,~,~,explained,~] = pca(resp_all2d');


    top_comp = score(:,1:num_comp);
    top_comp2 = reshape(top_comp, [num_t, num_tn, num_comp]);
    
    tomp_comp_all{n_gr} = top_comp2;
    exp_var_all{n_gr} = explained;
end


if num_comp >= 3
    figure; hold on;
    for n_tn = 1:num_tn
        tn1 = tn_all(n_tn);
        if sum(tn1 == [28, 19, 30])
            symb1 = '*';
        else
            symb1 = 'o';
        end
        plot3(top_comp2(:,n_tn, 1), top_comp2(:,n_tn, 2), top_comp2(:,n_tn, 3), 'color', app.ops.context_types_all_colors2{tn1}, 'LineWidth', 2);
        plot3(top_comp2(1,n_tn, 1), top_comp2(1,n_tn, 2), top_comp2(1,n_tn, 3), symb1, 'color', app.ops.context_types_all_colors2{tn1}, 'LineWidth', 2);
    end
    title(sprintf('comp 1-3; %s, region %s', title_tag, reg_tag));
end

if num_comp >= 6
    figure; hold on;
    for n_tn = 1:num_tn
        tn1 = tn_all(n_tn);
        if sum(tn1 == [28, 19, 30])
            symb1 = '*';
        else
            symb1 = 'o';
        end
        plot3(top_comp2(:,n_tn, 4), top_comp2(:,n_tn, 5), top_comp2(:,n_tn, 6), 'color', app.ops.context_types_all_colors2{tn1}, 'LineWidth', 2);
        plot3(top_comp2(1,n_tn, 4), top_comp2(1,n_tn, 5), top_comp2(1,n_tn, 6), symb1, 'color', app.ops.context_types_all_colors2{tn1}, 'LineWidth', 2);
    end
    title(sprintf('comp 4-6; %s, region %s', title_tag, reg_tag));
end


dist_all = cell(num_dist,1);
has_dist_idx = false(num_dist,1);
for n_list = 1:num_dist
    if sum(sum(tn_all == dist_list(n_list,:)')) == 2

        has_dist_idx(n_list) = 1;
        A = squeeze(top_comp2(:,tn_all == dist_list(n_list,1),:));
        B = squeeze(top_comp2(:,tn_all == dist_list(n_list,2),:));

        %dist1 = diag(pdist2(A,B,'euclidean'))
        dist1 = sum((A - B).^2,2).^(1/2);
        dist_all{n_list} = dist1;
    end
end
figure; hold on;
for n_list = 1:num_dist
    if has_dist_idx(n_list)
        plot(dist_all{n_list});
    end
end
title(sprintf('distance D to C; %s, region %s, %dcomp; %.2f%%var', title_tag, reg_tag, num_comp, sum(explained(1:num_comp))));
legend(dist_lab(has_dist_idx))

if plot_extra
    figure; 
    plot(explained(1:50))
    title(sprintf('explained variance  %s, region %s', title_tag, reg_tag));
    for n_comp = 1:num_comp
        figure; hold on
        for n_tn = 1:num_tn
            plot(squeeze(top_comp2(:,n_tn, n_comp)), 'color', app.ops.context_types_all_colors2{tn_all(n_tn)})
        end
        title(sprintf('comp %d', n_comp))
    end
end
end