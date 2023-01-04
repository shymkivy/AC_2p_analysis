function f_dv_plot_trial_ave_pca(app)

normalize1 = 0;
plot_extra = app.plotsuperdeetsCheckBox.Value;

num_comp = 6;

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
    
    reg_cell_labels = f_dv_get_area_label(app, data1);
    
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

num_effdsets = num_dsets*num_tn_gr;

resp_all2 = reshape(resp_all, num_effdsets, num_reg);
cell_counts2  = reshape(cell_counts, num_effdsets, num_reg);
mouse_id2 = reshape(mouse_id, num_effdsets, 1);
dset_id2 = reshape(dset_id, num_effdsets, 1);

[gr_id, gr_tag] = f_dv_combine_data(app, mouse_id2, dset_id2);


%% now the pca

gr_all = unique(gr_id);
num_gr = numel(gr_all);

top_comp_comb = cell(1, num_reg);
exp_var_comb = cell(1, num_reg);
for n_reg = 1:num_reg
    num_cells = sum(cell_counts2(:,n_reg));
    
    resp2 = cat(1,resp_all2{:,n_reg});
    
    resp2d = reshape(resp2, num_cells, []);
    
    if normalize1
        resp2d = f_normalize(resp2d, 'norm_mean_std');
    end

    %[U,S,V] = svd(resp_all2dn);
    % score*coeff'
    [~,score,~,~,explained,~] = pca(resp2d');
    
    top_comp = score(:,1:num_comp);
    top_comp2 = reshape(top_comp, [num_t, num_tn, num_comp]);
    
    top_comp_comb{1, n_reg} = top_comp2;
    exp_var_comb{1, n_reg} = explained;
end

colors_tn = app.ops.context_types_all_colors2;
tn1 = tn_all(1, :);

for n_reg = 1:num_reg
    top_comp2 = top_comp_comb{1, n_reg};
    exp_var2 = exp_var_comb{1, n_reg};
    if ~isempty(top_comp2)
        if num_comp >= 3
            trs1 = [1 2 3];
            sum_var = sum(exp_var2(trs1(1):trs1(3)));
            title_tag2 = sprintf('%s; combined; region %s; comp %d-%d; %.2f%% var', title_tag, leg_list{n_reg}, trs1(1), trs1(3), sum_var);
            f_dv_plot3_pc3(top_comp2(:,:,trs1), tn1, title_tag2, plot_t, colors_tn);
        end

        if 0%num_comp >= 6
            trs1 = [4 5 6];
            sum_var = sum(exp_var2(trs1(1):trs1(3)));
            title_tag2 = sprintf('%s; combined; region %s; comp %d-%d; %.2f%% var', title_tag, leg_list{n_reg}, trs1(1), trs1(3), sum_var);
            f_dv_plot3_pc3(top_comp2(:,:,trs1), tn1, title_tag2, plot_t, colors_tn);
        end
    end
end

%%
top_comp_all = cell(num_gr, num_reg);
exp_var_all = cell(num_gr, num_reg);
num_cells_all = zeros(num_gr, num_reg);

for n_gr = 1:num_gr
    idx_dset = gr_all(n_gr) == gr_id;
    for n_reg = 1:num_reg

        num_cells = sum(cell_counts2(idx_dset,n_reg));
        num_cells_all(n_gr, n_reg) = num_cells;
        
        if num_cells > 10
            resp2 = cat(1,resp_all2{idx_dset,n_reg});

            resp2d = reshape(resp2, num_cells, []);

            if normalize1
                resp2d = f_normalize(resp2d, 'norm_mean_std');
            end

            %[U,S,V] = svd(resp_all2dn);
            % score*coeff'
            [~,score,~,~,explained,~] = pca(resp2d');

            top_comp = score(:,1:num_comp);
            top_comp2 = reshape(top_comp, [num_t, num_tn, num_comp]);

            top_comp_all{n_gr, n_reg} = top_comp2;
            exp_var_all{n_gr, n_reg} = explained;
        end
    end
end

if plot_extra
    for n_gr = 1:num_gr
        for n_reg = 1:num_reg
            top_comp2 = top_comp_all{n_gr, n_reg};
            exp_var2 = exp_var_all{n_gr, n_reg};
            if ~isempty(top_comp2)
                if num_comp >= 3
                    trs1 = [1 2 3];
                    sum_var = sum(exp_var2(trs1(1):trs1(3)));
                    title_tag2 = sprintf('%s; %s n= %d; region %s; comp %d-%d; %.2f%% var', title_tag, gr_tag, n_gr, leg_list{n_reg}, trs1(1), trs1(3), sum_var);
                    f_dv_plot3_pc3(top_comp2(:,:,trs1), tn1, title_tag2, plot_t, colors_tn);
                end

                if 0%num_comp >= 6
                    trs1 = [4 5 6];
                    sum_var = sum(exp_var2(trs1(1):trs1(3)));
                    title_tag2 = sprintf('%s; %s n= %d; region %s; comp %d-%d; %.2f%% var', title_tag, gr_tag, n_gr, leg_list{n_reg}, trs1(1), trs1(3), sum_var);
                    f_dv_plot3_pc3(top_comp2(:,:,trs1), tn1, title_tag2, plot_t, colors_tn);
                end
            end
        end
    end
end

dist_all2 = cell(num_dist, 1);
has_dist_idx = false(num_dist, 1);
for n_list = 1:num_dist
    if sum(sum(tn1 == dist_list(n_list,:)')) == 2
        has_dist_idx(n_list) = 1;
        dist_all = cell(num_gr, num_reg);
        for n_gr = 1:num_gr
            for n_reg = 1:num_reg
                top_comp2 = top_comp_all{n_gr, n_reg};
                if ~isempty(top_comp2)
                    
                    A = squeeze(top_comp2(:,tn1 == dist_list(n_list,1),:));
                    B = squeeze(top_comp2(:,tn1 == dist_list(n_list,2),:));

                    %dist1 = diag(pdist2(A,B,'euclidean'))
                    dist1 = sum((A - B).^2,2).^(1/2);
                    dist_all{n_gr, n_reg} = dist1;
                end
            end
        end
        dist_all2{n_list} = dist_all;
    end
end

dist_all3 = dist_all2(has_dist_idx);
num_dist2 = numel(dist_all3);
dist_lab2 = dist_lab(has_dist_idx);


%%
onset_win = params.stats.onset_resp_win;
offset_win = params.stats.offset_resp_win;
mid_win = params.stats.middle_resp_win; %mid_win = [.3 .6];

win_frames{1} = logical((plot_t >= onset_win(1)) .* (plot_t <= onset_win(2)));
win_frames{2} = logical((plot_t >= offset_win(1)) .* (plot_t <= offset_win(2)));
win_frames{3} = logical((plot_t >= mid_win(1)) .* (plot_t <= mid_win(2))); 
labl1 = [app.ops.context_types_labels(tn_all(1,:))]; %, [{'Cont comb'; 'Red comb pool'; 'Dev comb'}]];
win_labels = {'Onset', 'Offset', 'Middle'};
num_win = numel(win_frames);


dist_mean_all = cell(num_dist2, num_reg);
dist_sem_all = cell(num_dist2, num_reg);

for n_list = 1:num_dist2
    dist_all4 = dist_all3{n_list};
    for n_reg = 1:num_reg
        dist_all5 = cat(2, dist_all4{:, n_reg});
        num_groups = size(dist_all5,2);
        
        dist_mean_all{n_list, n_reg} = mean(dist_all5,2);
        dist_sem_all{n_list, n_reg} = std(dist_all5, [], 2)/sqrt(max(num_groups-1, 1));
    end
end

dist_mean_allcat = cat(1, dist_mean_all{:});
dist_sem_allcat = cat(1, dist_sem_all{:});

ylim1 = ([min([dist_mean_allcat - dist_sem_allcat; 0]), max(dist_mean_allcat + dist_sem_allcat)*1.1]);



for n_list = 1:num_dist2
    dist_all4 = dist_all3{n_list};
    win_vals = cell(num_win,num_reg);
    reg_lab = cell(1,num_reg);
    
    figure; hold on; axis tight
    sall = cell(num_reg,1);
    for n_reg = 1:num_reg
        
        dist_all5 = cat(2, dist_all4{:, n_reg});
        
        for n_win = 1:num_win
            win_vals{n_win, n_reg} = mean(dist_all5(win_frames{n_win},:),1)';
        end
        
        num_groups = size(dist_all5,2);
        reg_lab{n_reg} = repmat(n_reg,num_groups,1);
        
        dist_mean = mean(dist_all5,2);
        dist_sem = std(dist_all5, [], 2)/sqrt(max(num_groups-1, 1));
        
        color1 = f_dv_get_leg_color(app, leg_list{n_reg});
        s1 = shadedErrorBar_YS(plot_t, dist_mean, dist_sem, color1);
        sall{n_reg} = s1.mainLine;
        ylim(ylim1);
        %plot(plot_t, dist_all5, 'color', color1, 'LineWidth', 2);
        
    end
    ylabel('euclidean distance');
    xlabel('Time')
    title(sprintf('distance %s; %s, region %s, %dcomp; %.2f%%var; ', dist_lab2{n_list}, title_tag, reg_tag, num_comp, sum(explained(1:num_comp))));
    legend([sall{:}], leg_list);
    
    reg_lab2 = cat(1, reg_lab{:});
    for n_win = 1:num_win
        win_data = cat(1, win_vals{n_win, :});
        
        [p1, tbl1, stats2]  = anova1(win_data, reg_lab2, 'off');
        title_tag4 = sprintf('%s; %s; %s reg; %s win', title_tag, dist_lab2{n_list}, reg_tag, win_labels{n_win});
        f_dv_plot_anova1(p1, tbl1, stats2, title_tag4);
        
    end
    
end


%%
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