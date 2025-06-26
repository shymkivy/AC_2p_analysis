function f_dv_plot_trial_ave_pca(app)

ops = app.ops;
params = f_dv_gather_params(app);
data0 = app.data;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);
[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);

normalize1 = 0;

num_comp = params.plot_pca_dim;
num_pl_d = str2double(app.plotaxesDropDown.Value);
num_pl = floor(num_comp/num_pl_d);
num_plot_comp = num_pl*num_pl_d;

trs_all = reshape(1:num_plot_comp, num_pl_d, []);

cell_num_lim = 1;

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
num_dsets = numel(data.experiment);

[plot_t, trial_frames] = f_dv_compute_window_t(params.trial_window, data.proc_data{1}.frame_data.volume_period_ave);

num_t = sum(trial_frames);

tn0 = f_dv_get_trial_number(params, data.MMN_freq{1});
tn00 = tn0(1,:);

[num_tn_gr, num_tn] = size(tn0);

%reg_all = app.ops.regions_to_analyze;
[region_num_all, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);

num_reg = size(region_num_all,1);

mouse_id = cell(num_tn_gr, num_dsets);
dset_id = zeros(num_tn_gr, num_dsets);

resp_all = cell(num_tn_gr, num_dsets, num_reg);
cell_counts = zeros(num_tn_gr, num_dsets, num_reg);

fprintf('dset #/%d: ', num_dsets);
for n_dset = 1:num_dsets
    fprintf('..%d', n_dset);
    data1 =  data(n_dset,:);
    stats1 = data1.stats{params.planes};
    params.n_dset = find(data1.idx == data0);

    cdata = f_dv_compute_cdata(data1, params);
    
    tn1 = f_dv_get_trial_number(params, data1.MMN_freq{1});

    %num_cells = sum([stats1.num_cells]);

    firing_rate = cdata.S_sm;
    trial_types = data1.trial_types{1};
    stim_times = data1.stim_frame_index{params.n_pl};
    mmn_freq = data1.MMN_freq{1};
    
    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
    [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, ops);
    
    reg_cell_labels = f_dv_get_area_label(data1, params, ops);
    
    for n_tngr = 1:num_tn_gr
        mouse_id{n_tngr, n_dset} = data1.mouse_id{1};
        dset_id(n_tngr, n_dset) = n_dset;

        resp_cells = f_dv_get_resp_vals_cells(stats1, tn1(n_tngr,:), params);
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
                    temp_resp = trial_data_sort_wctx(resp_cell_idx,:,(trial_types_wctx == ops.context_types_all(tn1(n_tngr, n_tn))));
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

[gr_id, gr_tag] = f_dv_combine_data(mouse_id2, dset_id2, params);


%% now the pca

gr_all = unique(gr_id);
num_gr = numel(gr_all);

top_comp_comb = cell(1, num_reg);
exp_var_comb = cell(1, num_reg);
pr_all = zeros(1, num_reg);
for n_reg = 1:num_reg
    num_cells = sum(cell_counts2(:,n_reg));
    if num_cells
        resp2 = cat(1,resp_all2{:,n_reg});
        
        resp2d = reshape(resp2, num_cells, []);
        
        if normalize1
            resp2d = f_normalize(resp2d, 'norm_mean_std');
        end
    
        %[U,S,V] = svd(resp_all2dn);
        % score*coeff'
        %warning('off', 'stats:pca:ColRankDefX')
        [~,score,~,~,explained,~] = pca(resp2d');
        
        top_comp = score(:,1:num_comp);
        top_comp2 = reshape(top_comp, [num_t, num_tn, num_comp]);
        
        top_comp_comb{1, n_reg} = top_comp2;
        exp_var_comb{1, n_reg} = explained;
        pr_all(1,n_reg) = (sum(exp_var_comb{1, n_reg}).^2)/sum(exp_var_comb{1, n_reg}.^2);
    end
end


colors_tn = ops.context_types_all_colors2;
%tn1 = tn_all(1, :);

for n_reg = 1:num_reg
    top_comp2 = top_comp_comb{1, n_reg};
    exp_var2 = exp_var_comb{1, n_reg};
    if ~isempty(top_comp2)
        for n_pl = 1:num_pl
            trs1 = trs_all(:, n_pl);
            if num_comp >= trs1(end)
                sum_var = sum(exp_var2(trs1(1):trs1(end)));
                title_tag2 = sprintf('%s; combined; region %s; comp %d-%d; %.2f%% var; pr=%.2f', title_tag, leg_list{n_reg}, trs1(1), trs1(end), sum_var, pr_all(1,n_reg));
                if num_pl_d == 3
                    f_dv_plot3_t_pc(top_comp2(:,:,trs1), tn00, title_tag2, colors_tn, params.shadow_on3d, params.shadow_axis_locs, params.render_painters, params.grid_on, params.reverse_xyz);
                elseif num_pl_d == 2
                    figure(); hold on;
                    for n_tn = 1:num_tn
                        plot(top_comp2(:,n_tn,trs1(1)), top_comp2(:,n_tn,trs1(2)), 'color', colors_tn{tn00(n_tn)});
                    end
                    xlabel(sprintf('PC%d', trs1(1)));
                    ylabel(sprintf('PC%d', trs1(2)));
                    title(title_tag2, 'interpreter', 'none');
                end
            end
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

if params.plot_super_deets
    for n_gr = 1:num_gr
        for n_reg = 1:num_reg
            top_comp2 = top_comp_all{n_gr, n_reg};
            exp_var2 = exp_var_all{n_gr, n_reg};
            if ~isempty(top_comp2)
                if num_comp >= 3
                    trs1 = [1 2 3];
                    sum_var = sum(exp_var2(trs1(1):trs1(3)));
                    title_tag2 = sprintf('%s; %s n= %d; region %s; comp %d-%d; %.2f%% var', title_tag, gr_tag, n_gr, leg_list{n_reg}, trs1(1), trs1(3), sum_var);
                    f_dv_plot3_t_pc(top_comp2(:,:,trs1), tn00, title_tag2, colors_tn, params.shadow_on3d, params.shadow_axis_locs, params.render_painters, params.grid_on, params.reverse_xyz);
                end

                if 0%num_comp >= 6
                    trs1 = [4 5 6];
                    sum_var = sum(exp_var2(trs1(1):trs1(3)));
                    title_tag2 = sprintf('%s; %s n= %d; region %s; comp %d-%d; %.2f%% var', title_tag, gr_tag, n_gr, leg_list{n_reg}, trs1(1), trs1(3), sum_var);
                    f_dv_plot3_t_pc(top_comp2(:,:,trs1), tn00, title_tag2, colors_tn, params.shadow_on3d, params.shadow_axis_locs, params.render_painters, params.grid_on, params.reverse_xyz);
                end
            end
        end
    end
end
%%
if params.pca_distances
    dist_all2 = cell(num_dist, 1);
    has_dist_idx = false(num_dist, 1);
    for n_list = 1:num_dist
        if sum(sum(tn00 == dist_list(n_list,:)')) == 2
            has_dist_idx(n_list) = 1;
            dist_all = cell(num_gr, num_reg);
            for n_gr = 1:num_gr
                for n_reg = 1:num_reg
                    top_comp2 = top_comp_all{n_gr, n_reg};
                    if ~isempty(top_comp2)
                        
                        A = squeeze(top_comp2(:,tn00 == dist_list(n_list,1),:));
                        B = squeeze(top_comp2(:,tn00 == dist_list(n_list,2),:));
    
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
    if strcmpi(params.stim_window, 'Onset')
        win1 = params.stats.resp_win_onset;
    elseif strcmpi(params.stim_window, 'Offset')
        win1 = params.stats.resp_win_offset;
    elseif strcmpi(params.stim_window, 'Middle')
        win1 = params.stats.resp_win_middle;
    end
    win_frames = logical((plot_t >= win1(1)) .* (plot_t <= win1(2)));

    
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
    
    labl2 = ops.context_types_labels(tn1);
    colors1 = ops.context_types_all_colors2(tn1);

    for n_list = 1:num_dist2
        dist_all4 = dist_all3{n_list};
        win_vals = cell(1,num_reg);
        reg_lab = cell(1,num_reg);
        
        figure; hold on; axis tight
        if params.plot_stim
            r1 = rectangle('Position', [0 ylim1(1) 0.5 ylim1(2)-ylim1(1)]);
            r1.FaceColor = [ops.context_types_all_colors2{params.stim_freq_color} params.stim_transparancy];
            r1.EdgeColor = [ops.context_types_all_colors2{params.stim_freq_color} params.stim_transparancy];
        end

        sall = cell(num_reg,1);
        for n_reg = 1:num_reg
            if num_cells_all(n_reg) >= cell_num_lim
                dist_all5 = cat(2, dist_all4{:, n_reg});
                
                win_vals{n_reg} = mean(dist_all5(win_frames,:),1)';
                
                num_groups = size(dist_all5,2);
                reg_lab{n_reg} = repmat(n_reg,num_groups,1);
                
                

                dist_mean = mean(dist_all5,2);
                dist_sem = std(dist_all5, [], 2)/sqrt(max(num_groups-1, 1));
                
                color1 = f_dv_get_leg_color(ops, leg_list{n_reg});
                s1 = shadedErrorBar_YS(plot_t, dist_mean, dist_sem, color1);
                sall{n_reg} = s1.mainLine;
                ylim(ylim1);
                %plot(plot_t, dist_all5, 'color', color1, 'LineWidth', 2);
            end
        end
        ylabel('euclidean distance');
        xlabel('Time')
        title(sprintf('distance %s; %s, region %s, %dcomp; %.2f%%var; ', dist_lab2{n_list}, title_tag, reg_tag, num_comp, sum(explained(1:num_comp))), 'Interpreter','none');
        legend([sall{:}], leg_list(logical(sum(num_cells_all,1))));
        
        title_tag4 = sprintf('%s; %s; %s reg; %s win', title_tag, dist_lab2{n_list}, reg_tag, win_label);

        f_plot_bar(win_vals, leg_list', ops.cond_colors, colors1)
        title(title_tag4, 'interpreter', 'none');
        ylim(ylim1);

        reg_lab2 = cat(1, reg_lab{:});
        win_data = cat(1, win_vals{:});
        
        [p1, tbl1, stats2]  = anova1(win_data, reg_lab2, 'off');
        
        f_dv_plot_anova1(p1, tbl1, stats2, title_tag4, leg_list);
        
    end
end
%%
max_plot_comp = 20;
figure; hold on; axis tight
for n_reg = 1:num_reg
    if num_cells_all(n_reg)
        max_plot_comp2 = min(max_plot_comp, num_cells_all(n_reg));
        color1 = f_dv_get_leg_color(ops, leg_list{n_reg});
        plot(0:max_plot_comp2 , [100; 100 - cumsum(exp_var_all{n_reg}(1:max_plot_comp2))], 'color', color1, 'LineWidth', 2);
    end
end
l1 = line([num_comp, num_comp], [0 100]);
l1.Color = [.5 .5 .5];
l1.LineStyle = '--';
ylabel('Residual variance');
xlabel('num components');
legend([leg_list(logical(sum(num_cells_all,1))), {'num comp used'}]);
title(sprintf('Residual variance  %s, region %s', title_tag, reg_tag), 'Interpreter','none');

if params.plot_super_deets
    figure; hold on;
    for n_reg = 1:num_reg
        if num_cells_all(n_reg)
            color1 = f_dv_get_leg_color(ops, leg_list{n_reg});
            plot(exp_var_all{n_reg}(1:50), 'color', color1);
        end
    end
    ylabel('explained variance');
    legend(leg_list(logical(sum(num_cells_all,1))));
    title(sprintf('explained variance  %s, region %s', title_tag, reg_tag), 'Interpreter','none');

end

end