function f_dv_plot_trial_ave_pca(app)

normalize1 = 0;
plot_extra = 0;

num_comp = 10;

dist_list = [18, 20;
             28, 30;
             19, 20;
             29, 30;
             18, 19;
             28, 29;
             2, 3;
             3, 4;
             4, 5;
             5, 6;
             6, 7;
             7, 8;
             8, 9];
         
dist_lab = {'DD-cont', 'DD-cont flip', 'DD-red', 'DD-red flip', 'Cont-red', 'Cont-red flip',...
            '2-3', '3-4', '4-5', '5-6', '6-7', '7-8', '8-9'};

num_dist = size(dist_list,1);

%%
n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = numel(data.experiment);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

num_t = sum(trial_frames);

tn_all = f_dv_get_trial_number(app);
[num_tn_gr, num_tn] = size(tn_all);

params = f_dv_gather_params(app);

reg_all = app.ops.regions_to_analyze;
[region_num_all, reg_tag, leg_list] = f_dv_get_region_sel_val2(app);

num_reg = size(region_num_all,1);

mouse_id = cell(num_tn_gr, num_dsets);
dset_id = zeros(num_tn_gr, num_dsets);
resp_all = cell(num_tn_gr, num_dsets, num_reg);
cell_counts = zeros(num_tn_gr, num_dsets, num_reg);

fprintf('dset #/%d: ', num_dsets);
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
    
    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
    [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);
    
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
    
    for n_tngr = 1:num_tn_gr
        mouse_id{n_tngr, n_dset} = data1.mouse_id{1};
        dset_id(n_tngr, n_dset) = n_dset;
        tn1 = tn_all(n_tngr,:);
        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn1);
        %cell_is_resp = stats1.peak_resp_cells(:,tn_all);

        for n_reg = 1:num_reg
            region_num = region_num_all(n_reg,:);

            reg_cell_idx = logical(sum(reg_cell_labels == region_num,2));

            % get resp cells
            resp_cell_idx = logical(sum(resp_cells,2).*reg_cell_idx);
            num_cells = sum(resp_cell_idx);

            if num_cells
                cell_counts(n_tngr, n_dset, n_reg) = num_cells;
                resp_all{n_tngr, n_dset, n_reg} = zeros(num_cells, num_t, num_tn);
                for n_tn = 1:num_tn
                    temp_resp = trial_data_sort_wctx(resp_cell_idx,:,(trial_types_wctx == app.ops.context_types_all(tn1(n_tn))));
                    resp_all{n_tngr, n_dset, n_reg}(:,:,n_tn) = mean(temp_resp,3);
                end
            end
        end
    end
end
fprintf('\n');

resp_all2 = reshape(resp_all, num_dsets*num_tn_gr, num_reg);
cell_counts2  = reshape(cell_counts, num_dsets*num_tn_gr, num_reg);
mouse_id2 = reshape(mouse_id, num_dsets*num_tn_gr, 1);
dset_id2 = reshape(dset_id, num_dsets*num_tn_gr, 1);
%% now the pca

if strcmpi(app.statsbetweenDropDown.Value, 'combine')
    
elseif strcmpi(app.statsbetweenDropDown.Value, 'Mice')
    groups1 = mouse_id;
elseif strcmpi(app.statsbetweenDropDown.Value, 'Dsets')
end

mice_all = unique(mouse_id);
num_mice = numel(mice_all);

top_comp_all = cell(num_mice, num_reg);
exp_var_all = cell(num_mice, num_reg);
num_cells_all = zeros(num_mice, num_reg);

for n_ms = 1:num_mice
    idx_dset = strcmpi(mice_all{n_ms}, mouse_id);
    for n_reg = 1:num_reg

        num_cells = sum(cell_counts(idx_dset,n_reg));
        num_cells_all(n_ms, n_reg) = num_cells;
        
        if num_cells > 10
            resp2 = cat(1,resp_all{idx_dset,n_reg});

            resp2d = reshape(resp2, num_cells, []);

            if normalize1
                resp2d = f_normalize(resp2d, 'norm_mean_std');
            end

            %[U,S,V] = svd(resp_all2dn);
            % score*coeff'
            [~,score,~,~,explained,~] = pca(resp2d');

            top_comp = score(:,1:num_comp);
            top_comp2 = reshape(top_comp, [num_t, num_tn, num_comp]);

            top_comp_all{n_ms, n_reg} = top_comp2;
            exp_var_all{n_ms, n_reg} = explained;
        end
    end
end

colors_tn = app.ops.context_types_all_colors2;
tn1 = tn_all(:,1);

for n_ms = 1:num_mice
    for n_reg = 1:num_reg
        top_comp2 = top_comp_all{n_ms, n_reg};
        exp_var2 = exp_var_all{n_ms, n_reg};
        if ~isempty(top_comp2)

            if num_comp >= 3
                trs1 = [1 2 3];
                sum_var = sum(exp_var2(trs1(1):trs1(3)));
                title_tag2 = sprintf('%s; mouse %s; region %s; %.2f var', title_tag, mice_all{n_ms}, leg_list{n_reg}, sum_var);
                f_dv_plot3_pc3(top_comp2, tn1, trs1, title_tag2, plot_t, colors_tn);
            end

            if 0%num_comp >= 6
                trs1 = [4 5 6];
                sum_var = sum(exp_var2(trs1(1):trs1(3)));
                title_tag2 = sprintf('%s; mouse %s; region %s; %.2f var', title_tag, mice_all{n_ms}, leg_list{n_reg}, sum_var);
                f_dv_plot3_pc3(top_comp2, tn1, trs1, title_tag2, plot_t, colors_tn);
            end
        end
    end
end

dist_all2 = cell(num_dist, 1);
has_dist_idx = false(num_dist, 1);
for n_list = 1:num_dist
    if sum(sum(tn1 == dist_list(n_list,:)')) == 2
        has_dist_idx(n_list) = 1;
        dist_all = cell(num_mice, num_reg);
        for n_ms = 1:num_mice
            for n_reg = 1:num_reg
                top_comp2 = top_comp_all{n_ms, n_reg};
                if ~isempty(top_comp2)
                    
                    A = squeeze(top_comp2(:,tn1 == dist_list(n_list,1),:));
                    B = squeeze(top_comp2(:,tn1 == dist_list(n_list,2),:));

                    %dist1 = diag(pdist2(A,B,'euclidean'))
                    dist1 = sum((A - B).^2,2).^(1/2);
                    dist_all{n_ms, n_reg} = dist1;
                end
            end
        end
        dist_all2{n_list} = dist_all;
    end
end

dist_all3 = dist_all2(has_dist_idx);
num_dist2 = numel(dist_all3);

for n_list = 1:num_dist2
    dist_all4 = dist_all3{n_list};
    figure; hold on; axis tight
    sall = cell(num_reg,1);
    for n_reg = 1:num_reg
        
        dist_all5 = cat(2, dist_all4{:, n_reg});
        num_groups = size(dist_all5,2);
        
        dist_mean = mean(dist_all5,2);
        dist_sem = std(dist_all5, [], 2)/sqrt(max(num_groups-1, 1));
        
        color1 = f_dv_get_leg_color(app, leg_list{n_reg});
        s1 = shadedErrorBar_YS(plot_t, dist_mean, dist_sem, color1);
        sall{n_reg} = s1.mainLine;
        %plot(plot_t, dist_all5, 'color', color1, 'LineWidth', 2);
    end
    ylabel('euclidean distance');
    xlabel('Time')
    title(sprintf('distance %s; %s, region %s, %dcomp; %.2f%%var; ', dist_lab{n_list}, title_tag, reg_tag, num_comp, sum(explained(1:num_comp))));
    legend([sall{:}], leg_list)
end

max_plot_comp = 20;
figure; hold on;
for n_reg = 1:num_reg
    color1 = f_dv_get_leg_color(app, leg_list{n_reg});
    plot(0:max_plot_comp , [100; 100 - cumsum(exp_var_all{n_reg}(1:max_plot_comp))], 'color', color1, 'LineWidth', 2)
end
l1 = line([num_comp, num_comp], [0 100]);
l1.Color = [.5 .5 .5];
l1.LineStyle = '--';
ylabel('Residual variance');
xlabel('num components');
legend([leg_list, {'num comp used'}]);
title(sprintf('Residual variance  %s, region %s', title_tag, reg_tag));

if plot_extra
    figure; hold on;
    for n_reg = 1:num_reg
        color1 = f_dv_get_leg_color(app, leg_list{n_reg});
        plot(exp_var_all{n_reg}(1:50), 'color', color1);
    end
    ylabel('explained variance');
    legend(leg_list);
    title(sprintf('explained variance  %s, region %s', title_tag, reg_tag));

    for n_comp = 1:num_comp
        figure; hold on
        for n_tn = 1:num_tn
            plot(squeeze(top_comp2(:,n_tn, n_comp)), 'color', app.ops.context_types_all_colors2{tn1(n_tn)})
        end
        title(sprintf('comp %d', n_comp))
    end
end
end