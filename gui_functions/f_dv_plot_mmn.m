function f_dv_plot_mmn(app)

z_ylim_max = 10;
z_ylim_min = -0.5;

add_combined = 1;

n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = numel(data.experiment);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

num_t = sum(trial_frames);

ctx_plot_list = [18, 19, 20; ...
                     28, 29, 30]';
                 
num_flip = size(ctx_plot_list,2);
                 
params = f_dv_gather_params(app);

[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val2(app);
num_plots = size(region_num,1);

resp_all = cell(num_dsets, num_flip, num_plots);
for n_reg_pl = 1:num_plots
    region_num2 = region_num(n_reg_pl,:);
    cell_counts = zeros(num_dsets, num_flip);
    for n_flip = 1:num_flip
        fprintf('Red %d/%d; Flip %d/%d; dset #/%d: ', n_reg_pl, num_plots, n_flip, num_flip, num_dsets);
        tn_all = ctx_plot_list(:,n_flip);
        ctx1 = app.ops.context_types_all(tn_all)';

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

            resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all);
            %cell_is_resp = stats1.peak_resp_cells(:,tn_all);

            reg_cell_labels = f_dv_get_area_label(app, data1);

            reg_cell_idx = logical(sum(reg_cell_labels == region_num2,2));

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

            % get resp cells
            resp_cell_idx = logical(sum(resp_cells,2).*reg_cell_idx);
            num_cells = sum(resp_cell_idx);
            cell_counts(n_dset, n_flip) = num_cells;

            resp_all{n_dset, n_flip, n_reg_pl} = zeros(num_cells, num_t, size(ctx1,2));
            for n_ctx = 1:size(ctx1,2)
                ctx2 = ctx1(:,n_ctx)';
                temp_resp = trial_data_sort_wctx(resp_cell_idx,:,logical(sum(trial_types_wctx == ctx2,2)));
                resp_all{n_dset, n_flip, n_reg_pl}(:,:,n_ctx) = mean(temp_resp,3);
                if sum(sum(isnan(resp_all{n_dset, n_flip}(:,:,n_ctx))))
                    1;
                end
            end
        end
        fprintf('\n');
    end
    
end
resp_all_pool = cell(num_plots, num_flip);
for n_reg_pl = 1:num_plots
    for n_flip = 1:num_flip
        resp_all_pool{n_reg_pl, n_flip} = cat(1,resp_all{:,n_flip, n_reg_pl});
    end
end

if add_combined
    for n_reg_pl = 1:num_plots
        resp_all_pool{n_reg_pl, num_flip+1} = cat(1,resp_all_pool{n_reg_pl, :});
    end
end

num_flip2 = size(resp_all_pool,2);

num_cells = zeros(num_plots, num_flip2);
resp_mean = cell(num_plots, num_flip2);
resp_sem = cell(num_plots, num_flip2);
y_lim_max = 0;
y_lim_min = 0;
for n_reg_pl = 1:num_plots
    for n_flip = 1:num_flip2
        temp_resp = resp_all_pool{n_reg_pl, n_flip};
        temp_num_cells = size(temp_resp,1);
        num_cells(n_reg_pl, n_flip) = temp_num_cells;
        resp_mean{n_reg_pl, n_flip} = squeeze(mean(temp_resp,1));
        resp_sem{n_reg_pl, n_flip} = squeeze(std(temp_resp,[],1)/sqrt(max(temp_num_cells-1,1)));

        max_vals = resp_mean{n_reg_pl, n_flip} + resp_sem{n_reg_pl, n_flip};
        min_vals = resp_mean{n_reg_pl, n_flip} - resp_sem{n_reg_pl, n_flip};
        y_lim_max = max([y_lim_max max(max_vals(:))]);
        y_lim_min = min([y_lim_min min(min_vals(:))]);
    end
end
y_lim_max = max([y_lim_max z_ylim_max]);
y_lim_min = min([y_lim_min z_ylim_min]);

for n_reg_pl = 1:num_plots
    figure;
    for n_flip = 1:num_flip2
        subplot(num_flip2,1,n_flip); hold on; axis tight; ylim([y_lim_min, y_lim_max]);
        if num_cells(n_reg_pl, n_flip)
            for n_ctx = 1:size(ctx_plot_list,1)
                color2 = app.ops.context_types_all_colors2{ctx_plot_list(n_ctx,1)};
                shadedErrorBar_YS(plot_t, resp_mean{n_reg_pl, n_flip}(:,n_ctx),resp_sem{n_reg_pl, n_flip}(:,n_ctx), color2);
            end

        end
        if app.ConverttoZCheckBox.Value
            ylabel('Z-score');
        else
            ylabel('response mag');
        end
        xlabel('Time (sec)');
        title(sprintf('%s; %d cells', reg_tag, num_cells(n_flip)));
    end
    sgtitle([title_tag '; region ' leg_list{n_reg_pl}], 'Interpreter', 'none');
end


%% plot differences

onset_win = params.stats.onset_resp_win;
offset_win = params.stats.offset_resp_win;
mid_win = [.3 .6];
win_frames{1} = logical((plot_t >= onset_win(1)) .* (plot_t <= onset_win(2)));
win_frames{2} = logical((plot_t >= offset_win(1)) .* (plot_t <= offset_win(2)));
win_frames{3} = logical((plot_t >= mid_win(1)) .* (plot_t <= mid_win(2))); 
labl1 = [app.ops.context_types_labels(ctx_plot_list), [{'Cont comb'; 'Red comb pool'; 'Dev comb'}]];
win_labels = {'Onset', 'Offset', 'Middle'};

% cont red dev
ylim1 = [-1 9];
for n_win = 1:3
    figure;
    for n_flip = 1:num_flip2
        means_all = zeros(num_plots, 3);
        sem_all = zeros(num_plots, 3);
        for n_reg_pl = 1:num_plots
            temp_mean = resp_mean{n_reg_pl, n_flip};
            temp_sem = resp_sem{n_reg_pl, n_flip};

            means_all(n_reg_pl,:) = mean(temp_mean(win_frames{n_win},:),1);
            sem_all(n_reg_pl, :) = mean(temp_sem(win_frames{n_win},:),1);
        end

        subplot(1, num_flip2, n_flip); hold on;
        bar(categorical(labl1(:,n_flip), labl1(:,n_flip)), [0 0 0]);
        for n1 = 1:3
            b1 = bar(n1, means_all(:,n1));
            errorbar(n1 + ((1:4)-2.5)/5.5, means_all(:,n1), sem_all(:,n1), '.k');
            for n_br = 1:numel(b1)
                b1(n_br).FaceColor = app.ops.cond_colors{n_br};
            end
        end
        ylim(ylim1)
        title(sprintf('%s resp flip %d', win_labels{n_win}, n_flip));
        ylabel('z-score');
    end
end


%% also differences

% cont red dev
categories1 = {'DD - red', 'DD - cont', 'Cont - red'};
sub_pairs = [3, 2; 3, 1; 1, 2];
for n_win = 1:3
    figure;
    for n_flip = 1:num_flip2
        means_all = zeros(num_plots, 3);
        sem_all = zeros(num_plots, 3);
        for n_reg_pl = 1:num_plots
            temp_mean = resp_mean{n_reg_pl, n_flip};
            temp_sem = resp_sem{n_reg_pl, n_flip};

            for n_pair = 1:3
                means_all(n_reg_pl,n_pair) = mean(temp_mean(win_frames{n_win},sub_pairs(n_pair,1)),1) - mean(temp_mean(win_frames{n_win},sub_pairs(n_pair,2)),1);
                sem_all(n_reg_pl,n_pair) = mean([mean(temp_sem(win_frames{n_win},sub_pairs(n_pair,1)),1), mean(temp_sem(win_frames{n_win},sub_pairs(n_pair,2)),1)]);
            end
        end

        subplot(1, num_flip2, n_flip); hold on;
        bar(categorical(categories1, categories1), [0 0 0]);
        for n1 = 1:3
            b1 = bar(n1, means_all(:,n1));
            errorbar(n1 + ((1:4)-2.5)/5.5, means_all(:,n1), sem_all(:,n1), '.k');
            for n_br = 1:numel(b1)
                b1(n_br).FaceColor = app.ops.cond_colors{n_br};
            end
        end
        ylim(ylim1)
        title(sprintf('%s resp flip %d', win_labels{n_win}, n_flip));
        ylabel('z-score');
    end
end


disp('Done');
end