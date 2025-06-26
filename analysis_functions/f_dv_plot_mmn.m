function f_dv_plot_mmn(data0, params, ops)

[data, title_tag] = f_dv_get_data_by_mouse_selection(data0, params);
num_dsets = numel(data.experiment);

[plot_t, trial_frames] = f_dv_compute_window_t(params.trial_window, data.proc_data{1}.frame_data.volume_period_ave);
num_t = sum(trial_frames);

%tn_all = [18, 19, 20; 28, 29, 30]';
if sum(strcmpi(params.trial_type, {'Context', 'Context_flip'}))
    leg1 = ops.context_types_labels;
else
    params.trial_type = 'Context_both_comb';
    leg1 = ops.context_types_labels_trim2;
end
tn_all = f_dv_get_trial_number(params);
[num_gr, num_tn] = size(tn_all);

[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);
num_reg = size(region_num,1);

mouse_id = cell(num_gr, num_dsets);
dset_id = zeros(num_gr, num_dsets);

resp_all = cell(num_dsets, num_gr, num_reg);
cell_counts = zeros(num_dsets, num_gr, num_reg);
fprintf('dset #/%d: ', num_dsets);
for n_dset = 1:num_dsets
    
    fprintf('..%d', n_dset);
    data1 =  data(n_dset,:);
    params.n_dset = find(data1.idx == data0.idx);
    
    [cdata, stats1] = f_dv_get_new_cdata_stats(data1, params);
    num_cells = sum([cdata.num_cells]);
    firing_rate = cat(1,cdata.S_sm);
    mmn_freq = data1.MMN_freq{1};
    trial_types = data1.trial_types{1};
    trial_types_ctx = f_dv_mark_tt_ctx(trial_types, mmn_freq, ops);

    stim_times = data1.stim_frame_index{params.n_pl};
    
    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
    %[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

    if params.convert_to_z
        st_mean_mean = cat(1,stats1.stat_trials_mean_mean);
        st_mean_sem = cat(1,stats1.stat_trials_mean_sem);
    else
        st_mean_mean = zeros(num_cells,1);
        st_mean_sem = ones(num_cells,1);
    end
    trial_data_sort = (trial_data_sort - st_mean_mean)./st_mean_sem;
    
    reg_cell_labels = f_dv_get_area_label(data1, params, ops);
    
    for n_gr = 1:num_gr
        mouse_id{n_gr, n_dset} = data1.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;
        
        tn1 = tn_all(n_gr,:);
        
        resp_cells = f_dv_get_resp_vals_cells(stats1, tn1, params);
        resp_cells2 = sum(resp_cells,2);
        %cell_is_resp = stats1.peak_resp_cells(:,tn_all);

        for n_reg = 1:num_reg
            n_reg2 = region_num(n_reg,:);
            
            % get resp cells
            reg_cell_idx = logical(sum(reg_cell_labels == n_reg2,2));
            resp_cell_idx = logical(resp_cells2.*reg_cell_idx);
            
            num_cells = sum(resp_cell_idx);

            if num_cells
                cell_counts(n_dset, n_gr, n_reg) = num_cells;
                resp_all{n_dset, n_gr, n_reg} = zeros(num_cells, num_t, num_tn);
                for n_ctx = 1:num_tn
                    ctx1 = ops.context_types_all(tn1(n_ctx));
                    %idx1 = trial_types_wctx == ctx2;
                    idx1 = logical(sum([trial_types, trial_types_ctx] == ctx1,2));
                    temp_resp = trial_data_sort(resp_cell_idx,:,idx1);
                    resp_all{n_dset, n_gr, n_reg}(:,:,n_ctx) = mean(temp_resp,3);
                end
            end
        end
    end
end
fprintf('\n');

num_effdsets = num_dsets*num_gr;

resp_all2 = reshape(resp_all, num_effdsets, num_reg);
cell_counts2  = reshape(cell_counts, num_effdsets, num_reg);
mouse_id2 = reshape(mouse_id, num_effdsets, 1);
dset_id2 = reshape(dset_id, num_effdsets, 1);

[gr_id, gr_tag] = f_dv_combine_data(mouse_id2, dset_id2, params);

gr_all = unique(gr_id);
num_gr2 = numel(gr_all);

%%
z_ylim_max = 0;
z_ylim_min = 0;

y_lim_max = 0;
y_lim_min = 0;

win_label = params.stim_window;

title_tag3 = sprintf('%s; %s; reg %s; %s', title_tag, params.responsive_cells_select, reg_tag, params.trial_type);
if params.convert_to_z
    title_tag3 = sprintf('%s; Z', title_tag3);
end
title_tag4 = sprintf('%s win; %s', win_label, title_tag3);

if strcmpi(win_label, 'Onset')
    win1 = params.stats.resp_win_onset;
elseif strcmpi(win_label, 'Offset')
    win1 = params.stats.resp_win_offset;
elseif strcmpi(win_label, 'Middle')
    win1 = params.stats.resp_win_middle;
end
win_frames = logical((plot_t >= win1(1)) .* (plot_t <= win1(2)));

resp_all3 = cell(num_reg, 1);
resp_mean = cell(num_reg,1);
resp_sem = cell(num_reg,1);
num_cells = zeros(num_reg,1);
resp_reg_tn_all = cell(num_reg, num_tn);
reg_lab = cell(num_reg,1);
for n_reg = 1:num_reg
    temp_resp = cat(1,resp_all2{:, n_reg});
    resp_all3{n_reg} = temp_resp;
    temp_num_cells = size(temp_resp,1);
    if temp_num_cells
        num_cells(n_reg, 1) = temp_num_cells;
        temp_mean = reshape(mean(temp_resp,1), num_t, num_tn);
        temp_sem = reshape(std(temp_resp,[],1)/sqrt(max(temp_num_cells-1,1)), num_t, num_tn);
        resp_mean{n_reg, 1} = temp_mean;
        resp_sem{n_reg, 1} = temp_sem;
        
        reg_lab{n_reg} = repmat(n_reg, temp_num_cells, 1);
        for n_tn = 1:num_tn
            resp_reg_tn_all{n_reg, n_tn} = mean(temp_resp(:,win_frames,n_tn),2);
        end
    
        max_vals = temp_mean + temp_sem;
        min_vals = temp_mean - temp_sem;
        y_lim_max = max([y_lim_max max(max_vals(:))]);
        y_lim_min = min([y_lim_min min(min_vals(:))]);
    end
end

y_lim_max = max([max([y_lim_max z_ylim_max])*1.05, params.max_y_lim]);
y_lim_min = min([min([y_lim_min z_ylim_min])*1.05, params.min_y_lim]);
ylim1 = [y_lim_min, y_lim_max];

tn1 = tn_all(1,:);
labl1 = ops.context_types_labels_trim2(tn1);

for n_reg = 1:num_reg
    leg_all1 = cell(num_tn,1);
    
    if num_cells(n_reg)
        figure; hold on; axis tight; ylim(ylim1);
        if params.plot_stim
            r1 = rectangle('Position', [0 y_lim_min 0.5 y_lim_max-y_lim_min]);
            r1.FaceColor = [ops.context_types_all_colors2{params.stim_freq_color} params.stim_transparancy];
            r1.EdgeColor = [ops.context_types_all_colors2{params.stim_freq_color} params.stim_transparancy];
        end
        for n_ctx = 1:num_tn
            color2 = ops.context_types_all_colors2{tn1(n_ctx)};
            sh1 = shadedErrorBar_YS(plot_t, resp_mean{n_reg, 1}(:,n_ctx),resp_sem{n_reg, 1}(:,n_ctx), color2);
            leg_all1{n_ctx} = sh1.mainLine;
        end
        if params.convert_to_z
            ylabel('Z-score');
        else
            ylabel('response mag');
        end
        xlabel('Time (sec)');
        legend([leg_all1{:}], leg1(tn1))
        title(sprintf('%s; %s; %d cells', leg_list{n_reg}, title_tag3, num_cells(n_reg)), 'Interpreter', 'none');
    end
end

%% plot differences
if params.plot_super_deets
    % cont red dev
    %ylim1 = [-1 9];

    resp_all4 = cell(num_reg, num_tn);
    for n_reg = 1:num_reg
        for n_tn = 1:num_tn
            resp_all4{n_reg, n_tn} = mean(resp_all3{n_reg}(:,win_frames,n_tn),2);
        end
    end
    colors1 = ops.context_types_all_colors2(tn1);
    labl2 = ops.context_types_labels(tn1);

    if strcmpi(params.bar_plot_type, 'bar')
        f_plot_bar(resp_all4, labl2, colors1, ops.cond_colors)
        ylim(ylim1)
        title(title_tag4, 'interpreter', 'none');
        ylabel('z-score');
    elseif strcmpi(params.bar_plot_type, 'bar_points')
        f_plot_bar_points(resp_all4, labl2, ops.context_colors, ops.cond_colors, 'std', 0)
        title(title_tag4, 'interpreter', 'none');
        ylabel('z-score');
    elseif strcmpi(params.bar_plot_type, 'violin')
        f_plot_bar_points(resp_all4, labl2, ops.context_colors, ops.cond_colors, 'std', 1)
        title(title_tag4, 'interpreter', 'none');
        ylabel('z-score');
    end

    if params.plot_stats
        if sum(strcmpi(reg_tag, {'all', 'primary vs secondary'}))
            % between regions
            for n_tn = 1:num_tn
                [p_all, tbl_all, stats_all]  = anova1(cat(1,resp_reg_tn_all{:,n_tn}),cat(1,reg_lab{:}), 'off');
                title_tag5 = sprintf('between reg; %s; %s', title_tag4, labl1{n_tn});
                f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag5, leg_list, 1);
            end
        else
            resp_temp = cat(2,resp_reg_tn_all{n_reg,:});
            % between trial types
            
            if strcmpi(params.bar_plot_type, 'parametric')
                [p_all, tbl_all, stats_all]  = anova1(resp_temp,[], 'off');
                title_tag5 = sprintf('between ctx; %s; %s', title_tag4, leg_list{n_reg}); 
                f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag5, labl1, 1);
                
                f_dv_plot_ttest_comb(resp_temp, title_tag4, labl1);
            elseif strcmpi(params.bar_plot_type, 'nonparametric')
                f_dv_plot_signedrank(resp_temp, title_tag4, labl1);
            end
        end
    end
    
    
    %% also differences
    % cont red dev
    categories1 = {'DD - red', 'DD - cont', 'Cont - red'};
    sub_pairs = [3, 2; 3, 1; 1, 2];

    means_all = zeros(num_reg, 3);
    sem_all = zeros(num_reg, 3);
    for n_reg = 1:num_reg
        if num_cells(n_reg)
            temp_mean = resp_mean{n_reg, 1};
            temp_sem = resp_sem{n_reg, 1};
    
            for n_pair = 1:3
                means_all(n_reg,n_pair) = mean(temp_mean(win_frames,sub_pairs(n_pair,1)),1) - mean(temp_mean(win_frames,sub_pairs(n_pair,2)),1);
                sem_all(n_reg,n_pair) = mean([mean(temp_sem(win_frames,sub_pairs(n_pair,1)),1), mean(temp_sem(win_frames,sub_pairs(n_pair,2)),1)]);
            end
        end
    end

    figure; hold on;
    bar(categorical(categories1, categories1), [0 0 0]);
    for n1 = 1:3
        num_reg2 = numel(means_all(:,n1));
        b1 = bar(n1, means_all(:,n1));
        if num_reg2 > 1
            errorbar(n1 + ((1:4)-2.5)/5.5, means_all(:,n1), sem_all(:,n1), '.k');
            for n_br = 1:numel(b1)
                b1(n_br).FaceColor = ops.cond_colors{n_br};
            end
        else
            errorbar(n1, means_all(:,n1), sem_all(:,n1), '.k');
            b1.FaceColor = [.4 .4 .4];
        end
    end
    ylim(ylim1)
    title(sprintf('ctx diff; %s win; %s', win_label, title_tag3), 'interpreter', 'none');
    ylabel('z-score');
    
    if params.plot_stats
        if sum(strcmpi(reg_tag, {'all', 'primary vs secondary'}))
            for n_pair = 1:3
                [p_all, tbl_all, stats_all]  = anova1(cat(1,resp_reg_tn_all{:,sub_pairs(n_pair,1)}) - cat(1,resp_reg_tn_all{:,sub_pairs(n_pair,2)}),cat(1,reg_lab{:}), 'off');
                title_tag5 = sprintf('between pair; %s; %s', title_tag4, categories1{n_pair});
                f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag5, leg_list);
            end
        end
    end
    
    %%
    for n_reg = 1:num_reg
        if num_cells(n_reg)
            temp_data = cat(1,resp_all2{:, n_reg});
            temp_data2 = reshape(temp_data, [], num_t*num_tn);
            
            norm_met = 2;
            if norm_met == 1
                temp_data3 = (temp_data2 - mean(temp_data2(:)));
                temp_data3 = temp_data3/std(temp_data3(:));
            elseif norm_met == 2
                temp_data3 = (temp_data2 - mean(temp_data2,2));
                temp_data3 = temp_data3./std(temp_data3,[],2);
            elseif norm_met == 3
                temp_data3 = temp_data2/std(temp_data2(:));
            elseif norm_met == 4
                temp_data3 = (temp_data2 - mean(temp_data2(:)));
                temp_data3 = temp_data3/norm(temp_data3(:));
            elseif norm_met == 5
                temp_data3 = temp_data2;
                temp_data3 = temp_data3/std(temp_data3(:));
            elseif norm_met == 0
                temp_data3 = temp_data2;
            end
            
            %temp_data3 = temp_data;
            hc_params.plot_dist_mat = 1;
            hc_params.distance_metric = 'cosine';
            hc_params.method = 'ward';
            clust_out = f_hcluster_wrap(temp_data3, hc_params);
        
            %[~, idx1] = sort(mean(temp_data(idx2,:,3),2), 'descend');
            %[~, idx2] = sort(mean(temp_data(:,:,1),2), 'ascend');
            
            % if app.InvertcmapCheckBox.Value
            %     temp_data = -temp_data;
            % else
            %     temp_data = temp_data;
            % end
            max_val = max(temp_data(:));
            min_val = min(temp_data(:));
            figure;
            %clim12 = ([min(temp_data2(:)), max(temp_data2(:))]);
            for n_ctx = 1:num_tn
                subplot(1,num_tn,n_ctx);
                imagesc(plot_t, [], temp_data(clust_out.dend_order,:,n_ctx));
                clim([min_val max_val]*0.6)
            end
            colormap(app.ColormapDropDown.Value)
            sgtitle(sprintf('%s; region %s; flip%d; nmet%d', title_tag3, leg_list{n_reg}, n_gr,norm_met), 'Interpreter', 'none');
        end
    end
end
disp('Done');
end