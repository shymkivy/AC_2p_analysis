function f_dv_plot_pop_mmn(app)
plot_first_tr = 1;
combine_reg = 1;
plot_split_mmn = 1;
num_red = 10;
alpha1 = app.StimtranspEditField.Value;
freq_col = app.stimcolorSpinner.Value;
color1 = app.ops.context_types_all_colors2{freq_col};
win_on = f_str_to_array(app.stats_OnsetRespwinEditField.Value);
win_off = f_str_to_array(app.stats_OffsetRespwinEditField.Value);

fit_type = 'best'; % exp, linear, best

win_use = [0.1, 0.9];

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = numel(data.dset_name);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);
num_t = sum(trial_frames);

if sum(strcmpi(app.trialtypeDropDown.Value, {'Context', 'Context_flip', 'Context_both_comb'}))
    tt_input = app.trialtypeDropDown.Value;
else
    tt_input = 'Context_both_comb';
end
tn_all = f_dv_get_trial_number(app, tt_input);
[num_gr, num_tn] = size(tn_all);


vol_per = app.ddata.proc_data{1}.frame_data.volume_period_ave;
plot_frames = round(num_red/vol_per*1000);
plot_t2 = ((1:plot_frames))*vol_per/1000;

params = f_dv_gather_params(app);

[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(app);
num_reg = size(region_num,1);

if num_reg>1
    col_reg = app.ops.cond_colors;
else
    col_reg = app.ops.context_colors(2);
end

mouse_id = cell(num_gr, num_dsets);
dset_id = zeros(num_gr, num_dsets);

num_cells_all = zeros(num_dsets, num_gr, num_reg);
data_first_trial = cell(num_dsets, num_gr, num_reg);
data_rem_trials_mean = cell(num_dsets, num_gr, num_reg);
data_rem_trials_cellmean = cell(num_dsets, num_gr, num_reg);
data_rem_trials_all = cell(num_dsets, num_gr, num_reg);
data_mmn_all = cell(num_dsets, num_gr, num_reg);
fprintf('dset #/%d: ', num_dsets);
for n_dset = 1:num_dsets
    fprintf('..%d', n_dset);
    data1 = data(n_dset,:);
    stats1 = data1.stats{n_pl};
    params.n_dset = find(data1.idx == app.data.idx);
    cdata = f_dv_compute_cdata(data1, params);
    
    if app.ConverttoZCheckBox.Value
        st_mean_mean = stats1.stat_trials_mean_mean;
        st_mean_sem = stats1.stat_trials_mean_sem;
    else
        st_mean_mean = zeros(cdata.num_cells,1);
        st_mean_sem = ones(cdata.num_cells,1);
    end

    firing_rate = cdata.S_sm;
    mmn_freq = data1.MMN_freq{1};
    trial_types = data1.trial_types{1};
    trial_types_ctx = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);
    stim_times = data1.stim_frame_index{n_pl};
    num_tt_all = numel(trial_types);

    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
    trial_data_sort = (trial_data_sort - st_mean_mean)./st_mean_sem;

    reg_cell_labels = f_dv_get_area_label(app, data1);

    for n_gr = 1:num_gr
        mouse_id{n_gr, n_dset} = data1.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;

        tn1 = tn_all(n_gr, :);
        num_tn = numel(tn1);
        red_start = app.ops.context_types_all(tn1(1) - 7);
        red_end = app.ops.context_types_all(50 - tn1(3));
        mmn_start = find(trial_types == red_start);
        num_tr = numel(mmn_start);
        mmn_end = zeros(num_tr,1);

        for n_st = 1:num_tr
            curr_ntt = mmn_start(n_st);
            tt1 = trial_types(curr_ntt);
            while and(tt1 >= red_start, tt1 < red_end)
                curr_ntt = curr_ntt + 1;
                if curr_ntt > num_tt_all
                    tt1 = red_end;
                else
                    tt1 = trial_types(curr_ntt);
                end
            end
            mmn_end(n_st) = curr_ntt-1;
        end

        frame_start = stim_times(mmn_start);
        frame_end = stim_times(mmn_end) + trial_frames(2);
        stim_dur = frame_end - frame_start+1;
        max_stim_dur = max([max(stim_dur), plot_frames]);

        trial_data_sort_red = nan(cdata.num_cells, max_stim_dur, num_tr);
        
        for n_tr = 1:num_tr
            trial_data_sort_red(:, 1:stim_dur(n_tr), n_tr) = firing_rate(:,frame_start(n_tr):frame_end(n_tr));
        end
        
        trial_data_sort_red = (trial_data_sort_red - st_mean_mean)./st_mean_sem;
        
        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn1);
        resp_cells2 = logical(sum(resp_cells,2));
        
        for n_reg = 1:num_reg
            n_reg2 = region_num(n_reg,:);
            
            % get resp cells
            reg_cell_idx = logical(sum(reg_cell_labels == n_reg2,2));
            resp_cell_idx = logical(resp_cells2.*reg_cell_idx);
            
            num_cells = sum(resp_cell_idx);
            num_cells_all(n_dset, n_gr, n_reg) = num_cells;
            
            if num_cells

                trial_data_sort_red2 = trial_data_sort_red(resp_cell_idx,:,:);
                trial_data_sort2 = trial_data_sort(resp_cell_idx,:,:);
                
                data_first_trial{n_dset, n_gr, n_reg} = trial_data_sort_red2(:, 1:plot_frames, 1);
                data_rem_trials_mean{n_dset, n_gr, n_reg} = mean(trial_data_sort_red2(:, 1:plot_frames, 2:end), 3, "omitnan");
                data_rem_trials_cellmean{n_dset, n_gr, n_reg} = mean(permute(trial_data_sort_red2(:, 1:plot_frames, 2:end), [3 2 1]), 3, "omitnan");
                data_rem_trials_all{n_dset, n_gr, n_reg} = trial_data_sort_red2(:, 1:plot_frames, 2:end);
                
                data_mmn_all{n_dset, n_gr, n_reg} = zeros(num_cells, num_t, num_tn);
                for n_tn = 1:num_tn
                    ctx1 = app.ops.context_types_all(tn1(n_tn));

                    idx1 = logical(sum([trial_types, trial_types_ctx] == ctx1,2));
                    temp_resp = trial_data_sort2(:, :, idx1);
                    data_mmn_all{n_dset, n_gr, n_reg}(:,:,n_tn) = mean(temp_resp,3);
                end
            end
        end
    end
end
fprintf('\n');

num_effdsets = num_dsets*num_gr;

data_first_trial2 = reshape(data_first_trial, num_effdsets, num_reg);
data_rem_trials_mean2 = reshape(data_rem_trials_mean, num_effdsets, num_reg);
data_rem_trials_cellmean2 = reshape(data_rem_trials_cellmean, num_effdsets, num_reg);
data_rem_trials_all2 = reshape(data_rem_trials_all, num_effdsets, num_reg);
data_mmn_all2 = reshape(data_mmn_all, num_effdsets, num_reg);

mouse_id2 = reshape(mouse_id, num_effdsets, 1);
dset_id2 = reshape(dset_id, num_effdsets, 1);

[gr_id, gr_tag] = f_dv_combine_data(app, mouse_id2, dset_id2);

gr_all = unique(gr_id);
num_gr2 = numel(gr_all);

%%



title_tag2 = sprintf('%s; %s; %s; %s reg', title_tag, app.ResponsivecellsselectDropDown.Value, tt_input, reg_tag);

n_gr = 1;
tn1 = tn_all(n_gr, :);
num_tn = numel(tn1);
ymin = 0;
ymax = 0;
% get ylims first
for n_reg = 1:num_reg
    first_tr_all = cat(1, data_first_trial2{:, n_reg});
    tr_mean_all = cat(1, data_rem_trials_mean2{:, n_reg});
    data_mmn_all3 = cat(1, data_mmn_all2{:, n_reg});
    
    num_cells = size(first_tr_all,1);
    if num_cells
        first_tr_mean  = mean(first_tr_all,1, "omitnan");
        first_tr_sem = std(first_tr_all, "omitnan")/sqrt(num_cells-1);
        
        tr_mean_mean = mean(tr_mean_all, 1, "omitnan");
        tr_mean_sem = std(tr_mean_all, "omitnan")/sqrt(num_cells-1);

        data_mmn_all3_mean = mean(data_mmn_all3, 1);
        data_mmn_all3_sem = std(data_mmn_all3, [], 1)/sqrt(num_cells-1);
        
        if plot_first_tr
            ymin = min([ymin, min([first_tr_mean-first_tr_sem, tr_mean_mean-tr_mean_sem, data_mmn_all3_mean(:)'-data_mmn_all3_sem(:)'])*1.05]);
            ymax = max([ymax, max([first_tr_mean+first_tr_sem, tr_mean_mean+tr_mean_sem, data_mmn_all3_mean(:)'+data_mmn_all3_sem(:)'])*1.05]);
        else
            ymin = min([ymin, min([tr_mean_mean-tr_mean_sem, data_mmn_all3_mean(:)'-data_mmn_all3_sem(:)'])*1.05]);
            ymax = max([ymax, max([tr_mean_mean+tr_mean_sem, data_mmn_all3_mean(:)'+data_mmn_all3_sem(:)'])*1.05]);
        end
    end
end
ymin = min([ymin, app.minYlimEditField.Value]);
ymax = max([ymax, app.maxYlimEditField.Value]);

if combine_reg
    if plot_first_tr
        f1 = figure(); axis tight; hold on
        if app.PlotstimCheckBox.Value
            if_plot_stim_rect(f1, num_red, [ymin, ymax], color1, alpha1)
        end
    end
    f2 = figure(); axis tight; hold on
    if app.PlotstimCheckBox.Value
        if_plot_stim_rect(f2, num_red, [ymin, ymax], color1, alpha1)
    end
end

pl_all1 = cell(num_reg,1);
pl_all2 = cell(num_reg,1);
for n_reg = 1:num_reg
    title_tag3 = sprintf('%s; %s', title_tag2, leg_list{n_reg});
    
    first_tr_all = cat(1, data_first_trial2{:, n_reg});
    tr_mean_all = cat(1, data_rem_trials_mean2{:, n_reg});

    data_mmn_all3 = cat(1, data_mmn_all2{:, n_reg});
    
    num_cells = size(first_tr_all,1);
    if num_cells
        
        first_tr_mean  = mean(first_tr_all,1, "omitnan");
        first_tr_sem = std(first_tr_all, "omitnan")/sqrt(num_cells-1);
        
        tr_mean_mean = mean(tr_mean_all, 1, "omitnan");
        tr_mean_sem = std(tr_mean_all, "omitnan")/sqrt(num_cells-1);

        data_mmn_all3_mean = mean(data_mmn_all3, 1);
        data_mmn_all3_sem = std(data_mmn_all3, [], 1)/sqrt(num_cells-1);
        
        if combine_reg
            if num_reg>1
                color2 = app.ops.cond_colors{n_reg};
            else
                color2 = app.ops.context_colors{2};
            end
        else
            color2 = app.ops.context_types_all_colors2{19};
        end
        
        if plot_first_tr
            if combine_reg
                figure(f1);
            else
                f1 = figure(); axis tight; hold on
                if app.PlotstimCheckBox.Value
                    if_plot_stim_rect(f1, num_red, [ymin, ymax], color1, alpha1)
                end
                title(sprintf('first trial %s; %d cells', title_tag3, num_cells), 'interpreter', 'none');
                if app.ConverttoZCheckBox.Value
                    ylabel('Z-score');
                else
                    ylabel('response mag');
                end
            end

            s1 = shadedErrorBar_YS(plot_t2, first_tr_mean, first_tr_sem, color2);
            pl_all1{n_reg} = s1.mainLine;
            ylim([ymin, ymax]);
            xlim([plot_t2(1), plot_t2(end)]);
            
        end
        if combine_reg
            figure(f2);
        else
            f2 = figure(); hold on;
            if app.PlotstimCheckBox.Value
                if_plot_stim_rect(f1, num_red, ylims1, color1, alpha1)
            end
            title(sprintf('red 1-%d trials %s; %d cells', num_red, title_tag3, num_cells), 'interpreter', 'none');
            if app.ConverttoZCheckBox.Value
                ylabel('Z-score');
            else
                ylabel('response mag');
            end
        end
        s1 = shadedErrorBar_YS(plot_t2, tr_mean_mean, tr_mean_sem, color2);
        pl_all2{n_reg} = s1.mainLine;
        ylim([ymin, ymax]);
        xlim([plot_t2(1), plot_t2(end)])
        
        if plot_split_mmn
            f3 = figure; 
            for n_tn = 1:num_tn
                tt1 = tn1(n_tn);
                subplot(1, 3, n_tn); hold on; axis tight;
                if_plot_stim_rect(f3, 1, [ymin, ymax], color1, alpha1)
                color2 = app.ops.context_types_all_colors2{tt1};
                shadedErrorBar_YS(plot_t, data_mmn_all3_mean(1,:,n_tn), data_mmn_all3_sem(1,:,n_tn), color2);
                ylim([ymin, ymax]);
                if app.ConverttoZCheckBox.Value
                    ylabel('Z-score');
                else
                    ylabel('response mag');
                end
            end
            sgtitle(sprintf('context trials %s; %d cells', title_tag3, num_cells), 'interpreter', 'none')
        end
    end
end

if combine_reg
    if plot_first_tr
        figure(f1);
        legend([pl_all1{:}], leg_list);
        title(sprintf('first trial %s', title_tag2), 'interpreter', 'none');
        if app.ConverttoZCheckBox.Value
            ylabel('Z-score');
        else
            ylabel('response mag');
        end
        xlabel('Time (sec)');
    end
    figure(f2);
    legend([pl_all2{:}], leg_list)
    title(sprintf('red 1-%d trials %s', num_red, title_tag2), 'interpreter', 'none');
    if app.ConverttoZCheckBox.Value
        ylabel('Z-score');
    else
        ylabel('response mag');
    end
    xlabel('Time (sec)');
end

%% do bar plot versions

if plot_first_tr
    f1 = figure(); hold on; axis tight;
    f1_fit_leg = cell(num_reg,1);
    f1_line = cell(num_reg, 1);
end
f2 = figure(); hold on; axis tight;
f2_fit_leg = cell(num_reg,1);
f2_line = cell(num_reg, 1);

first_tr_mean_all2 = zeros(num_reg, num_red);
first_tr_sem_all2 = zeros(num_reg, num_red);
tr_mean_all2 = zeros(num_reg, num_red);
tr_sem_all2 = zeros(num_reg, num_red);
tau_fit = zeros(num_reg, 1);
tau_fit_std = zeros(num_reg, 1);
tau_fit_df = zeros(num_reg, 1);
r_sq_all = zeros(num_reg, 1);
ft_tau_fit = zeros(num_reg, 1);
ft_tau_fit_std = zeros(num_reg, 1);
ft_tau_fit_df = zeros(num_reg, 1);
ft_r_sq_all = zeros(num_reg, 1);
for n_reg = 1:num_reg
    first_tr_all = cat(1, data_first_trial2{:, n_reg});
    tr_mean_all = cat(1, data_rem_trials_mean2{:, n_reg});

    num_cells = size(first_tr_all,1);
    if num_cells
        first_tr_mean_all3 = zeros(num_cells, num_red);
        tr_mean_all3 = zeros(num_cells, num_red);
        for n_red = 1:num_red
            win_idx = and(plot_t2 >= ((n_red-1)+win_use(1)), plot_t2 <= ((n_red-1)+win_use(2)));
            temp_mean1 = mean(first_tr_all(:,win_idx),2);
            
            first_tr_mean_all3(:,n_red) = temp_mean1;

            first_tr_mean_all2(n_reg, n_red) = mean(temp_mean1, 'omitnan');
            first_tr_sem_all2(n_reg, n_red) = std(temp_mean1, [], 1, 'omitnan')./sqrt(num_cells-1);
            
            temp_mean2 = mean(tr_mean_all(:,win_idx),2);

            tr_mean_all3(:,n_red) = temp_mean2;

            tr_mean_all2(n_reg, n_red) = mean(temp_mean2, 'omitnan');
            tr_sem_all2(n_reg, n_red) = std(temp_mean2, [], 1, 'omitnan')./sqrt(num_cells-1);
        end
        
        x_fit = 1:0.1:num_red;
        x_data = (1:num_red)';
        if plot_first_tr
            figure(f1);
            f1_line{n_reg} = errorbar(x_data, first_tr_mean_all2(n_reg,:), first_tr_sem_all2(n_reg,:), '.', color=col_reg{n_reg}, markersize=15);
            
            fit_data = if_get_fit_wrap(x_data, first_tr_mean_all2(n_reg,:)', fit_type);

            y_fit = fit_data.fit_eq(x_fit, fit_data.fit_pars);
            plot(x_fit, y_fit, '--', color=col_reg{n_reg});
                
            tau_mean = fit_data.fit_pars(2);
            tau_std = sqrt(fit_data.param_cov(2,2));
            ft_tau_fit(n_reg) = tau_mean;
            ft_tau_fit_std(n_reg) = tau_std;
            ft_tau_fit_df(n_reg) = fit_data.df;

            ft_r_sq_all(n_reg) = fit_data.r_sq_adj;
            if strcmpi(fit_data.fit_type, 'exp')
                leg_tag = sprintf('; tau=%.2f', fit_data.fit_pars(2));
            elseif strcmpi(fit_data.fit_type, 'linear')
                leg_tag = sprintf('; m=%.2f', fit_data.fit_pars(2));
            end

            f1_fit_leg{n_reg} = sprintf('%s; %s fit; R^2=%.2f%s', leg_list{n_reg}, fit_data.fit_type, fit_data.r_sq_adj, leg_tag);
            fprintf('%s first tr; fit eq: %s; param cov: %s\n', leg_list{n_reg}, fit_data.fit_eq_str, num2str(diag(fit_data.param_cov)'));
            fprintf('%s first tr; tau mean= %.3f; tau std=%.3f; df=%d\n', leg_list{n_reg}, tau_mean, tau_std, fit_data.df);
        end

        figure(f2);
        f2_line{n_reg} = errorbar(x_data, tr_mean_all2(n_reg,:), tr_sem_all2(n_reg,:), '.', color=col_reg{n_reg}, markersize=15);
        
        fit_data = if_get_fit_wrap(x_data, tr_mean_all2(n_reg,:)', fit_type);
        
        y_fit = fit_data.fit_eq(x_fit, fit_data.fit_pars);
        plot(x_fit, y_fit, '--', color=col_reg{n_reg});
        
        tau_mean = fit_data.fit_pars(2);
        tau_std = sqrt(fit_data.param_cov(2,2));
        tau_fit(n_reg) = tau_mean;
        tau_fit_std(n_reg) = tau_std;
        tau_fit_df(n_reg) = fit_data.df;
        r_sq_all(n_reg) = fit_data.r_sq_adj;
        
        if strcmpi(fit_data.fit_type, 'exp')
            leg_tag = sprintf('; tau=%.2f', fit_data.fit_pars(2));
        elseif strcmpi(fit_data.fit_type, 'linear')
            leg_tag = sprintf('; m=%.2f', fit_data.fit_pars(2));
        end

        f2_fit_leg{n_reg} = sprintf('%s; %s fit; R^2=%.2f%s', leg_list{n_reg}, fit_data.fit_type, fit_data.r_sq_adj, leg_tag);
        fprintf('%s trial ave; fit eq: %s; param cov: %s\n', leg_list{n_reg}, fit_data.fit_eq_str, num2str(diag(fit_data.param_cov)'));
        fprintf('%s trial ave; tau mean= %.3f; tau std=%.3f; df=%d\n', leg_list{n_reg}, tau_mean, tau_std, fit_data.df);
        if app.plotstatsCheckBox.Value
            % some anova, but a lot of groups, maybe better do some
            % regression
            % 
            % [p_all,tbl_all,stats_all] = anova1(first_tr_mean_all3(:),x1(:), 'nodisplay');
            % 
            % idx1 = strcmpi(tbl_all(1,:), 'MS');
            % idx2 = strcmpi(tbl_all(1,:), 'df');
            % idx5 = strcmpi(tbl_all(:,1), 'Error');
            % 
            % MSE = tbl_all{idx5, idx1};
            % dfe = tbl_all{idx5,idx2};
            % 
            % means1 = mean(first_tr_mean_all3, 1);
            % mean_tot = mean(first_tr_mean_all3(:));
            % 
            % t_vals1 = (means1 - mean_tot)/sqrt(MSE)*sqrt(num_cells-1);
            % p_vals1 = (1 - tcdf(abs(t_vals1), dfe))*2;
            % 
            % f_dv_plot_anova1(p_all, tbl_all, stats_all, leg_list{n_reg}, '');
            % some regresions
            

        end
    end
end
if plot_first_tr
    figure(f1);
    legend([f1_line{:}], f1_fit_leg);
    title(sprintf('first trials fits; %s', title_tag2), 'interpreter', 'none');
    if app.ConverttoZCheckBox.Value
        ylabel('Z-score');
    else
        ylabel('response mag');
    end
    xlabel('Time (sec)');
end

figure(f2);
legend([f2_line{:}], f2_fit_leg);
title(sprintf('red 1-%d fits; %s', num_red, title_tag2), 'interpreter', 'none');
if app.ConverttoZCheckBox.Value
    ylabel('Z-score');
else
    ylabel('response mag');
end
xlabel('Time (sec)');

if strcmpi(reg_tag, 'All comb')
    colors1 = {[.6 .6 .6]};
    colors2 = {'k'};
    colors3 = {'r'};
else
    colors1 = app.ops.cond_colors;
    colors2 = app.ops.cond_colors;
    colors3 = app.ops.cond_colors;
end

if plot_first_tr
    ft_tau_fit(ft_r_sq_all<0.1) = 0;
    ft_tau_fit_std2 = ft_tau_fit_std;
    ft_tau_fit_std2(ft_r_sq_all<0.1) = 0;

    figure; hold on;
    bar(categorical(leg_list,leg_list), ft_tau_fit, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
    for n_br = 1:num_reg
        br1 = bar(n_br, ft_tau_fit(n_br));
        br1.FaceColor = colors1{n_br};
    end
    if num_gr > 1
        errorbar(1:num_reg, ft_tau_fit, ft_tau_fit_std2, '.k','LineWidth',1);
    end
    title(sprintf('%s; first tr tau', title_tag2), 'Interpreter', 'none');
    ylabel('Tau (s)');
    
    f_dv_plot_param_ttest(ft_tau_fit, ft_tau_fit_std, ft_tau_fit_df, sprintf('%s; first tr tau', title_tag2), leg_list)
end

tau_fit(r_sq_all<0.1) = 0;
tau_fit_std2 = tau_fit_std;
tau_fit_std2(r_sq_all<0.1) = 0;
figure; hold on;
bar(categorical(leg_list,leg_list), tau_fit, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
for n_br = 1:num_reg
    br1 = bar(n_br, tau_fit(n_br));
    br1.FaceColor = colors1{n_br};
end
if num_gr > 1
    errorbar(1:num_reg, tau_fit, tau_fit_std2, '.k','LineWidth',1);
end
title(sprintf('%s; first tr tau', title_tag2), 'Interpreter', 'none');
ylabel('Tau (s)');

f_dv_plot_param_ttest(tau_fit, tau_fit_std, tau_fit_df, sprintf('%s; red 1-%d tau', title_tag2, num_red), leg_list)


% plot_bar_type = 2;
% if plot_bar_type == 1
%     if plot_first_tr
%         figure(); axis tight; hold on
%         b1 = bar(1:n_red, first_tr_mean_all2);
%         groupwidth = min(0.8, num_reg/(num_reg + 1.5));
%         for n_reg = 1:num_reg
%             b1(n_reg).FaceColor = app.ops.cond_colors{n_reg};
%             b1(n_reg).EdgeColor = app.ops.cond_colors{n_reg};
%             b1(n_reg).FaceAlpha = 0.4;
%             x = (1:num_red) - groupwidth/2 + (2*n_reg-1) * groupwidth / (2*num_reg);
%             errorbar(x, first_tr_mean_all2(n_reg,:), first_tr_sem_all2(n_reg,:), '.', color=app.ops.cond_colors{n_reg});
%         end
%         legend(b1, leg_list)
%     end
% 
%     figure(); axis tight; hold on
%     b1 = bar(1:n_red, tr_mean_all2);
%     groupwidth = min(0.8, num_reg/(num_reg + 1.5));
%     for n_reg = 1:num_reg
%         b1(n_reg).FaceColor = app.ops.cond_colors{n_reg};
%         b1(n_reg).EdgeColor = app.ops.cond_colors{n_reg};
%         b1(n_reg).FaceAlpha = 0.4;
%         x = (1:num_red) - groupwidth/2 + (2*n_reg-1) * groupwidth / (2*num_reg);
%         errorbar(x, tr_mean_all2(n_reg,:), tr_sem_all2(n_reg,:), '.', color=app.ops.cond_colors{n_reg});
%     end
%     legend(b1, leg_list)
% elseif plot_bar_type == 2
%     figure(); axis tight; hold on
%     for n_reg = 1:num_reg
%         errorbar(1:num_red, tr_mean_all2(n_reg,:), tr_sem_all2(n_reg,:), '.', color=app.ops.cond_colors{n_reg});
%     end
%     legend(leg_list)
% 
%     if plot_first_tr
%         figure(); axis tight; hold on
%         for n_reg = 1:num_reg
%             errorbar(1:num_red, first_tr_mean_all2(n_reg,:), first_tr_sem_all2(n_reg,:), '.', color=app.ops.cond_colors{n_reg});
%         end
%         legend(leg_list)
%     end
% end
%%
if app.plotsuperdeetsCheckBox.Value
    for n_reg = 1:num_reg
        title_tag3 = sprintf('%s; %s', title_tag2, leg_list{n_reg});
    
        data_rem_trials_cellmean3 = cat(1, data_rem_trials_cellmean2{:,n_reg});
    
        if ~isempty(data_rem_trials_cellmean3)    
            figure; 
            imagesc(data_rem_trials_cellmean3)
            title(sprintf('mean over cells; %s', title_tag3),  'interpreter', 'none')
        end
        data_rem_trials_mean3 = cat(1, data_rem_trials_mean2{:,n_reg});
        if ~isempty(data_rem_trials_mean3)
            figure; 
            imagesc(data_rem_trials_mean3)
            title(sprintf('mean over trials; %s', title_tag3),  'interpreter', 'none')
        end
    end
end

end

function  if_plot_stim_rect(fig, num_red, ylims1, color1, alpha1)

figure(fig);
for n_st = 1:num_red
    r1 = rectangle('Position', [n_st-1, ylims1(1), 0.5, ylims1(2)-ylims1(1)]);
    r1.FaceColor = [color1 alpha1];
    r1.EdgeColor = [color1 alpha1];
end

end

function fit_data = if_get_fit_wrap(data_x, data_y, fit_type)

if strcmpi(fit_type, 'exp')
    fit_data = f_get_fit(data_x, data_y, 'exp');
elseif strcmpi(fit_type, 'linear')
    fit_data = f_get_fit(data_x, data_y, 'linear');
elseif strcmpi(fit_type, 'best')
    fit_data_exp = f_get_fit(data_x, data_y, 'exp');
    fit_data_lin = f_get_fit(data_x, data_y, 'linear');
    total_fit = [fit_data_exp; fit_data_lin];

    [~, fit_idx] = max([total_fit.r_sq_adj]);

    fit_data = total_fit(fit_idx);
end

end
