function f_dv_plot_mmn(app)

z_ylim_max = 10;
z_ylim_min = -0.5;

%add_combined = 1;

plot_stats = app.plotstatsCheckBox.Value;
%%
n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = numel(data.experiment);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

num_t = sum(trial_frames);

%tn_all = [18, 19, 20; 28, 29, 30]';

if sum(strcmpi(app.trialtypeDropDown.Value, {'Context', 'Context_flip', 'Context_both_comb'}))
    tt_input = app.trialtypeDropDown.Value;
else
    tt_input = 'Context_both_comb';
end
tn_all = f_dv_get_trial_number(app, tt_input);
[num_gr, num_tn] = size(tn_all);

params = f_dv_gather_params(app);

onset_win = params.stats.onset_resp_win;
offset_win = params.stats.offset_resp_win;
mid_win = params.stats.middle_resp_win; %mid_win = [.3 .6];
win_frames{1} = logical((plot_t >= onset_win(1)) .* (plot_t <= onset_win(2)));
win_frames{2} = logical((plot_t >= offset_win(1)) .* (plot_t <= offset_win(2)));
win_frames{3} = logical((plot_t >= mid_win(1)) .* (plot_t <= mid_win(2))); 
win_labels = {'Onset', 'Offset', 'Middle'};
num_win = numel(win_frames);

[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val2(app);
num_reg = size(region_num,1);

mouse_id = cell(num_gr, num_dsets);
dset_id = zeros(num_gr, num_dsets);

resp_all = cell(num_dsets, num_gr, num_reg);
cell_counts = zeros(num_dsets, num_gr, num_reg);
fprintf('dset #/%d: ', num_dsets);
for n_dset = 1:num_dsets
    
    fprintf('..%d', n_dset);
    data1 =  data(n_dset,:);
    stats1 = data1.stats{n_pl};
    params.n_dset = find(data1.idx == app.data.idx);

    cdata = f_dv_compute_cdata(data1, params);

    firing_rate = cdata.S_sm;
    trial_types = data1.trial_types{1};
    stim_times = data1.stim_frame_index{n_pl};
    mmn_freq = data1.MMN_freq{1};

    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
    [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

    if app.ConverttoZCheckBox.Value
        st_mean_mean = stats1.stat_trials_mean_mean;
        st_mean_sem = stats1.stat_trials_mean_sem;
    else
        st_mean_mean = zeros(cdata.num_cells,1);
        st_mean_sem = ones(cdata.num_cells,1);
    end
    
    trial_data_sort_wctx = (trial_data_sort_wctx - st_mean_mean)./st_mean_sem;
    
    reg_cell_labels = f_dv_get_area_label(app, data1);
    
    for n_gr = 1:num_gr
        mouse_id{n_gr, n_dset} = data1.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;
        
        tn_all1 = tn_all(n_gr,:);
        ctx1 = app.ops.context_types_all(tn_all1)';
        
        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all1);
        %cell_is_resp = stats1.peak_resp_cells(:,tn_all);

        for n_reg = 1:num_reg
            n_reg2 = region_num(n_reg,:);
            
            % get resp cells
            reg_cell_idx = logical(sum(reg_cell_labels == n_reg2,2));
            resp_cell_idx = logical(sum(resp_cells,2).*reg_cell_idx);
            
            num_cells = sum(resp_cell_idx);
            cell_counts(n_dset, n_gr, n_reg) = num_cells;

            resp_all{n_dset, n_gr, n_reg} = zeros(num_cells, num_t, num_tn);
            for n_ctx = 1:num_tn
                ctx2 = ctx1(n_ctx);
                temp_resp = trial_data_sort_wctx(resp_cell_idx,:,trial_types_wctx == ctx2);
                resp_all{n_dset, n_gr, n_reg}(:,:,n_ctx) = mean(temp_resp,3);
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

[gr_id, gr_tag] = f_dv_combine_data(app, mouse_id2, dset_id2);

gr_all = unique(gr_id);
num_gr2 = numel(gr_all);

%%
y_lim_max = 0;
y_lim_min = 0;

resp_mean = cell(num_reg,1);
resp_sem = cell(num_reg,1);
num_cells = zeros(num_reg,1);
resp_reg_tn_all = cell(num_reg, num_tn, num_win);
reg_lab = cell(num_reg,1);
for n_reg = 1:num_reg
    
    temp_resp = cat(1,resp_all2{:, n_reg});
    temp_num_cells = size(temp_resp,1);
    num_cells(n_reg, 1) = temp_num_cells;
    temp_mean = reshape(mean(temp_resp,1), num_t, num_tn);
    temp_sem = reshape(std(temp_resp,[],1)/sqrt(max(temp_num_cells-1,1)), num_t, num_tn);
    resp_mean{n_reg, 1} = temp_mean;
    resp_sem{n_reg, 1} = temp_sem;
    
    reg_lab{n_reg} = repmat(n_reg, temp_num_cells, 1);
    for n_tn = 1:num_tn
        for n_win = 1:3
            resp_reg_tn_all{n_reg, n_tn, n_win} = mean(temp_resp(:,win_frames{n_win},n_tn),2);
        end
    end

    max_vals = temp_mean + temp_sem;
    min_vals = temp_mean - temp_sem;
    y_lim_max = max([y_lim_max max(max_vals(:))]);
    y_lim_min = min([y_lim_min min(min_vals(:))]);
end

y_lim_max = max([y_lim_max z_ylim_max]);
y_lim_min = min([y_lim_min z_ylim_min]);
ylim1 = [y_lim_min, y_lim_max];

title_tag3 = sprintf('%s; %s reg; %s', title_tag, reg_tag, tt_input);

tn_all1 = tn_all(1,:);
labl1 = app.ops.context_types_labels(tn_all1);

for n_reg = 1:num_reg
    leg_all1 = cell(num_tn,1);
    
    figure; hold on; axis tight; ylim(ylim1);
    if num_cells(n_reg, 1)
        for n_ctx = 1:num_tn
            color2 = app.ops.context_types_all_colors2{tn_all1(n_ctx)};
            sh1 = shadedErrorBar_YS(plot_t, resp_mean{n_reg, 1}(:,n_ctx),resp_sem{n_reg, 1}(:,n_ctx), color2);
            leg_all1{n_ctx} = sh1.mainLine;
        end

    end
    if app.ConverttoZCheckBox.Value
        ylabel('Z-score');
    else
        ylabel('response mag');
    end
    xlabel('Time (sec)');
    legend([leg_all1{:}], app.ops.context_types_labels(tn_all1))
    title(sprintf('%s; %s; %d cells', leg_list{n_reg}, title_tag3, num_cells(n_gr)), 'Interpreter', 'none');
    
    if plot_stats
        for n_win = 1:num_win
            [p_all, tbl_all, stats_all]  = anova1(cat(2,resp_reg_tn_all{n_reg,:,n_win}),[], 'off');
            title_tag4 = sprintf('between ctx; %s; %s win; %s', leg_list{n_reg}, win_labels{n_win}, title_tag3); 
            f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
        end
    end
end

%% plot differences

% cont red dev
%ylim1 = [-1 9];
for n_win = 1:3
    means_all = zeros(num_reg, 3);
    sem_all = zeros(num_reg, 3);
    
    for n_reg = 1:num_reg
        temp_mean = resp_mean{n_reg, 1};
        temp_sem = resp_sem{n_reg, 1};

        means_all(n_reg,:) = mean(temp_mean(win_frames{n_win},:),1);
        sem_all(n_reg, :) = mean(temp_sem(win_frames{n_win},:),1);
    end

    figure; hold on;
    bar(categorical(labl1(:,1), labl1(:,1)), [0 0 0]);
    for n1 = 1:3
        b1 = bar(n1, means_all(:,n1));
        errorbar(n1 + ((1:4)-2.5)/5.5, means_all(:,n1), sem_all(:,n1), '.k');
        for n_br = 1:numel(b1)
            b1(n_br).FaceColor = app.ops.cond_colors{n_br};
        end
    end
    ylim(ylim1)
    title(sprintf('%s resp flip %d', win_labels{n_win}, 1));
    ylabel('z-score');
    
    if plot_stats
        for n_tn = 1:num_tn
            [p_all, tbl_all, stats_all]  = anova1(cat(1,resp_reg_tn_all{:,n_tn,n_win}),cat(1,reg_lab{:}), 'off');
            title_tag4 = sprintf('between reg; %s win; %s; %s', win_labels{n_win}, labl1{n_tn}, title_tag3);
            f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
        end
    end
end


%% also differences

% cont red dev
categories1 = {'DD - red', 'DD - cont', 'Cont - red'};
sub_pairs = [3, 2; 3, 1; 1, 2];
for n_win = 1:3
    means_all = zeros(num_reg, 3);
    sem_all = zeros(num_reg, 3);
    for n_reg = 1:num_reg
        temp_mean = resp_mean{n_reg, 1};
        temp_sem = resp_sem{n_reg, 1};

        for n_pair = 1:3
            means_all(n_reg,n_pair) = mean(temp_mean(win_frames{n_win},sub_pairs(n_pair,1)),1) - mean(temp_mean(win_frames{n_win},sub_pairs(n_pair,2)),1);
            sem_all(n_reg,n_pair) = mean([mean(temp_sem(win_frames{n_win},sub_pairs(n_pair,1)),1), mean(temp_sem(win_frames{n_win},sub_pairs(n_pair,2)),1)]);
        end
    end

    figure; hold on;
    bar(categorical(categories1, categories1), [0 0 0]);
    for n1 = 1:3
        b1 = bar(n1, means_all(:,n1));
        errorbar(n1 + ((1:4)-2.5)/5.5, means_all(:,n1), sem_all(:,n1), '.k');
        for n_br = 1:numel(b1)
            b1(n_br).FaceColor = app.ops.cond_colors{n_br};
        end
    end
    ylim(ylim1)
    title(sprintf('%s resp flip %d', win_labels{n_win}, 1));
    ylabel('z-score');
    
    if plot_stats
        for n_pair = 1:3
            [p_all, tbl_all, stats_all]  = anova1(cat(1,resp_reg_tn_all{:,sub_pairs(n_pair,1),n_win}) - cat(1,resp_reg_tn_all{:,sub_pairs(n_pair,2),n_win}),cat(1,reg_lab{:}), 'off');
            title_tag4 = sprintf('between pair; %s win; %s; %s', win_labels{n_win}, categories1{n_pair}, title_tag3);
            f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
        end
    end
    
end

%%
for n_reg = 1:num_reg
 
    temp_data = cat(1,resp_all2{:, n_reg});

    hc_params.plot_dist_mat = 0;
    clust_out = f_hcluster_wrap(reshape(temp_data, [], num_t*num_tn), hc_params);

    %[~, idx1] = sort(mean(temp_data(idx2,:,3),2), 'descend');
    %[~, idx2] = sort(mean(temp_data(:,:,1),2), 'ascend');
    figure;
    clim1 = ([min(temp_data(:)), max(temp_data(:))]);
    for n_ctx = 1:num_tn
        subplot(1,num_tn,n_ctx);
        imagesc(plot_t, [], temp_data(clust_out.dend_order,:,n_ctx));
        caxis(clim1*0.6)
    end
    sgtitle(sprintf('%s; region %s; flip%d', title_tag3, leg_list{n_reg}, n_gr), 'Interpreter', 'none');
end

disp('Done');
end