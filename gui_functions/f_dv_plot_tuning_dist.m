function f_dv_plot_tuning_dist(app)

plot_ens = app.plotensamblesCheckBox.Value;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

if strcmpi(data(1,:).paradigm, 'tone_mmn')
    stim_params = data(end,:).proc_data{1}.stim_params;
    sf = stim_params.start_freq;
    incf = stim_params.increase_factor;

    ind2freq = @(ind) sf*(incf.^(ind-1));
else
    ind2freq = @(ind) ind;
    incf = 1/10;
end

tn_all = f_dv_get_trial_number(app);
[num_gr, num_tn] = size(tn_all);
num_dsets = numel(data.experiment);
reg_all = app.ops.regions_to_analyze;

[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(app);
num_reg = size(region_num,1);
leg_list2 = leg_list(region_num(:,1));

title_tag2 = sprintf('%s; %s reg', title_tag, reg_tag);

categories = app.ops.context_types_labels(tn_all(1,:));

stim_all = zeros(num_gr, num_dsets, num_tn);
for n_dset = 1:num_dsets
    data1 = data(n_dset,:);
    trial_types = data1.trial_types{1};
    trial_types_ctx = f_dv_mark_tt_ctx(trial_types, data1.MMN_freq{1}, app.ops);
    
    trial_types2 = [trial_types, trial_types_ctx];
    for n_gr = 1:num_gr
        stim_all(n_gr, n_dset, :) = squeeze(sum(sum(trial_types2 == reshape(app.ops.context_types_all(tn_all(n_gr,:)), 1, 1, []),1),2));
    end
end

stim_all2 = reshape(stim_all, num_dsets*num_gr, num_tn);
stim_all_mean = mean(stim_all2,1);
stim_all_sem = std(stim_all2, [], 1)./sqrt(max(num_dsets*num_gr - 1,1));
figure; hold on;
bar(categorical(categories,categories), stim_all_mean);
errorbar([], stim_all_mean, stim_all_sem, '.k')
title(sprintf('%s; stimuli counts', title_tag2), 'Interpreter', 'none');
ylabel('number of stimuli per dset');

%%
resp_cell_all = cell(num_gr, num_dsets);
num_cells_all = cell(num_gr, num_dsets);
loco_cells_all = cell(num_gr, num_dsets);
mouse_id = cell(num_gr, num_dsets);
dset_id = zeros(num_gr, num_dsets);

for n_dset = 1:num_dsets
    data1 = data(n_dset,:);

    if ~plot_ens
        stats1 = cat(1,data1.stats{n_pl});
    else
        stats1 = cat(1,data1.ensemble_tuning_stats{n_pl});
    end

    num_cells = sum([stats1.num_cells]);

    data_reg_avial = false;
    if app.UseregdatalabelsCheckBox.Value
        if ~isempty(data1.registered_data{1})
            reg_cell_labels = data1.registered_data{1}.reg_labels;
            data_reg_avial = true;
        else
            fprintf('dset %s reg data unavailable\n', data1.dset_name{1});
        end
    end
    if ~data_reg_avial
        reg_idx = find(strcmpi(reg_all, data1.area));
        reg_cell_labels = ones(num_cells,1)*reg_idx;
    end
    
    loco_cells = cat(1,[stats1.loco_resp_cells])';
    
        %reg_cell_idx = logical(sum(reg_cell_labels == region_num,2));
    for n_gr = 1:num_gr
        tn1 = tn_all(n_gr, :);
        mouse_id{n_gr, n_dset} = data1.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;
        
        [~, ~, ~, ~, resp_cells] = f_dv_get_resp_vals_cells(app, stats1, tn1);
        
        cell_is_resp_reg = zeros(num_cells, num_reg, num_tn);
        loco_cells_reg = zeros(num_cells, num_reg);
        num_cells1 = zeros(1, num_reg);
        for n_reg = 1:num_reg
            reg_idx = logical(sum(reg_cell_labels == region_num(n_reg,:),2));
            cell_is_resp_reg(reg_idx, n_reg, :) = resp_cells(reg_idx, :);
            loco_cells_reg(reg_idx, n_reg) = loco_cells(reg_idx);
            num_cells1(n_reg) = sum(reg_idx);
        end

        num_cells_all{n_gr, n_dset} = num_cells1;
        resp_cell_all{n_gr, n_dset} = cell_is_resp_reg;  
        
        loco_cells_all{n_gr, n_dset} = loco_cells_reg;
    end
end

num_effdsets = num_dsets*num_gr;
resp_cell_all2 = reshape(resp_cell_all, num_effdsets, 1);
loco_cells_all2 = reshape(loco_cells_all, num_effdsets, 1);
num_cells_all2  = reshape(num_cells_all, num_effdsets, 1);
mouse_id2 = reshape(mouse_id, num_effdsets, 1);
dset_id2 = reshape(dset_id, num_effdsets, 1);

[gr_id, gr_tag] = f_dv_combine_data(app, mouse_id2, dset_id2);

gr_id_uq = unique(gr_id);
num_gr = numel(gr_id_uq);
tuned_frac_all = cell(num_gr, 1);
tuned_frac_margtn_all = cell(num_gr, 1);
loco_frac_all = cell(num_gr, 1);
num_cells_gr = cell(num_gr, 1);

if num_tn == 10
    freq_lmh = cell(num_gr, 1);
end

for n_gr = 1:num_gr
    idx1 = gr_id_uq(n_gr) == gr_id;
    resp_cell_all3 = cat(1, resp_cell_all2{idx1});
    num_cells_all3 = sum(cat(1, num_cells_all2{idx1}),1);
    loco_cells_all3 = cat(1, loco_cells_all2{idx1});
    
    num_cells_gr{n_gr} = num_cells_all3;
    tuned_frac_all{n_gr} = sum(resp_cell_all3,1)./num_cells_all3;
    tuned_frac_margtn_all{n_gr} = sum(logical(sum(resp_cell_all3,3)),1)/sum(num_cells_all3);
    loco_frac_all{n_gr} = sum(loco_cells_all3,1)./num_cells_all3;
    if num_tn == 10
        resp3 = logical(cat(3,sum(resp_cell_all3(:,:,1:3),3), sum(resp_cell_all3(:,:,4:6),3), sum(resp_cell_all3(:,:,7:10),3)));
        freq_lmh{n_gr} = sum(resp3,1)./num_cells_all3;
    end
end

title_tag3 = sprintf('%s; stats by %s, %d groups', title_tag2, gr_tag, num_gr);

num_cells_gr2 = cat(1, num_cells_gr{:});
tuned_frac_all2 = cat(1, tuned_frac_all{:});
loco_frac_all2 = cat(1, loco_frac_all{:});
tuned_frac_margtn_all2 = cat(1, tuned_frac_margtn_all{:});
if num_tn == 10
    freq_lmh2 = cat(1,freq_lmh{:});
end

if strcmpi(reg_tag, 'All comb')
    colors1 = {[.6 .6 .6]};
    colors2 = {'k'};
    colors3 = {'r'};
else
    colors1 = app.ops.cond_colors;
    colors2 = app.ops.cond_colors;
    colors3 = app.ops.cond_colors;
end

min_cells = 5;

tun_mean = zeros(num_tn,num_reg);
tun_sem = zeros(num_tn,num_reg);
tun_count = zeros(1, num_reg);
tun_data_all = cell(num_reg,num_tn);
reg_lab_all = cell(num_reg, 1);
if num_tn == 10
    tun_lmh_mean = zeros(3, num_reg);
    tun_lmh_sem = zeros(3, num_reg);
    tun_lmh_all = cell(num_reg, 3);
end

for n_reg = 1:num_reg
    has_cells = num_cells_gr2(:,n_reg)>min_cells;
    num_gr_w_cells = sum(has_cells);
    tun_data2 = tuned_frac_all2(has_cells,n_reg,:);
    tun_mean(:,n_reg) = mean(tun_data2,1);
    tun_count(:,n_reg) = sum(has_cells);
    tun_sem(:,n_reg) = std(tun_data2,1)/sqrt(tun_count(:,n_reg) - 1);
    for n_tn = 1:num_tn
        tun_data_all{n_reg, n_tn} = tun_data2(:,1,n_tn);
    end
    reg_lab_all{n_reg} = repmat(n_reg, num_gr_w_cells, 1);
    if num_tn == 10
        tun_lmh_mean(:,n_reg) = mean(freq_lmh2(has_cells, n_reg,:),1);
        tun_lmh_sem(:,n_reg) = std(freq_lmh2(has_cells, n_reg,:),1)/sqrt(tun_count(:,n_reg) - 1);
        for n_tn = 1:3
            tun_lmh_all{n_reg, n_tn} = freq_lmh2(has_cells, n_reg,n_tn);
        end
    end
end

f1 = figure; hold on;
b1 = bar(categorical(categories,categories), tun_mean, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
for n_bar = 1:num_reg
    b1(n_bar).FaceColor = colors1{n_bar};
    b1(n_bar).EdgeColor = 'k';
    b1(n_bar).LineWidth = 1;
end
xb = nan(num_tn, num_reg);
if num_gr > 1
    for n_bar = 1:num_reg
        xb(:, n_bar) = b1(n_bar).XEndPoints;
    end
    errorbar(xb, tun_mean, tun_sem, '.k','LineWidth',1);
end
title(sprintf('%s; resp counts', title_tag3), 'Interpreter', 'none');
ylabel('Responsive cell fraction');
max_y = f1.Children.YLim;
f1.Children.YLim(2) = max([max_y, app.maxYlimEditField.Value]);
legend(b1, leg_list2)

if app.plotstatsCheckBox.Value
    for n_tn = 1:num_tn
        [p_tun, tbl_tun, stats_tun]  = anova1(cat(1, tun_data_all{:,n_tn}), cat(1, reg_lab_all{:}), 'off');
        title_tag4 = sprintf('resp counts; %s; %s', title_tag3, categories{n_tn});
        f_dv_plot_anova1(p_tun, tbl_tun, stats_tun, title_tag4, leg_list2);
    end
    if num_reg == 1
        tags = repmat(1:num_tn, [sum(has_cells),1 ]);
        [p_tun, tbl_tun, stats_tun]  = anova1(cat(1, tun_data_all{:}), tags(:), 'off');
        title_tag4 = sprintf('resp counts between trials; %s', title_tag3);
        f_dv_plot_anova1(p_tun, tbl_tun, stats_tun, title_tag4, categories,1);
    end
end

if num_tn == 10
    cat_lmh = {'freq 1-3', 'freq 4-6', 'freq 7-10'};
    f1 = figure; hold on;
    b1 = bar(categorical(cat_lmh,cat_lmh), tun_lmh_mean, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
    for n_bar = 1:num_reg
        b1(n_bar).FaceColor = colors1{n_bar};
        b1(n_bar).EdgeColor = 'k';
        b1(n_bar).LineWidth = 1;
    end
    xb = nan(3, num_reg);
    if num_gr > 1
        for n_bar = 1:num_reg
            xb(:, n_bar) = b1(n_bar).XEndPoints;
        end
        errorbar(xb, tun_lmh_mean, tun_lmh_sem, '.k','LineWidth',1);
    end
    title(sprintf('multi freq tuning; %s', title_tag3), 'Interpreter', 'none');
    ylabel('Responsive cell fraction');
    max_y = f1.Children.YLim;
    f1.Children.YLim(2) = max([max_y, app.maxYlimEditField.Value]);
    legend(b1, leg_list2)
    
    if app.plotstatsCheckBox.Value
        for n_tn = 1:3
            [p_tun, tbl_tun, stats_tun]  = anova1(cat(1, tun_lmh_all{:,n_tn}), cat(1, reg_lab_all{:}), 'off');
            title_tag4 = sprintf('multi freq tuning; %s; %s; tuning stats', title_tag3, cat_lmh{n_tn});
            f_dv_plot_anova1(p_tun, tbl_tun, stats_tun, title_tag4, leg_list2);
        end
    end

end

% do gauss fit
if app.plotstatsCheckBox.Value
    
    x_fit = (1:0.1:num_tn)';
    x_data = (1:num_tn)';
    
    f1 = figure(); hold on;
    f2 = figure(); hold on;
    line_f1 = cell(2, num_reg);
    line_f2 = cell(2, num_reg);
    leg_f1 = cell(2, num_reg);
    leg_f2 = cell(2, num_reg);
    for n_reg = 1:num_reg
        y_data = tun_mean(:,n_reg);
        y_sem = tun_sem(:,n_reg);

        fit_data_g = f_get_fit(x_data, y_data, 'gaus');
        fit_data_gl = f_get_fit(x_data, y_data, 'gaus+line', [fit_data_g.fit_pars; 0]);
        
        y_fit_g = fit_data_g.fit_eq(x_fit, fit_data_g.fit_pars);
        figure(f1)
        line_f1{1,n_reg} = plot(x_data, y_data, color=colors2{n_reg});
        errorbar(x_data, y_data, y_sem, '.','LineWidth',1, color=colors2{n_reg});
        line_f1{2,n_reg} = plot(x_fit, y_fit_g, '--', color=colors3{n_reg});
        leg_f1{1,n_reg} = sprintf('%s data', leg_list2{n_reg});
        leg_f1{2,n_reg} = sprintf('%s fit; R^2=%.2f;\\mu=%.2fkHz; \\sigma=%.2foct', leg_list2{n_reg}, fit_data_g.r_sq_adj, ind2freq(fit_data_g.fit_pars(2))/1000, fit_data_g.fit_pars(3)*incf);
        
        y_fit = fit_data_gl.fit_eq(x_fit, fit_data_gl.fit_pars);
        figure(f2);
        line_f2{1,n_reg} = plot(x_data, y_data, color=colors2{n_reg});
        errorbar(x_data, y_data, y_sem, '.','LineWidth',1, color=colors2{n_reg});
        line_f2{2,n_reg} = plot(x_fit, y_fit, '--', color=colors3{n_reg});
        leg_f2{1,n_reg} = sprintf('%s data', leg_list2{n_reg});
        leg_f2{2,n_reg} = sprintf('%s fit; R^2=%.2f;\\mu=%.2fkHz; \\sigma=%.2foct', leg_list2{n_reg}, fit_data_gl.r_sq_adj, ind2freq(fit_data_gl.fit_pars(2))/1000, fit_data_gl.fit_pars(3)*incf);
        
    end


    figure(f1);
    title(sprintf('gauss fit; %s', title_tag3), 'interpreter', 'none');
    legend([line_f1{:}], leg_f1(:))

    
    figure(f2);  
    title(sprintf('gauss+line fit; %s', title_tag3), 'interpreter', 'none');
    legend([line_f2{:}], leg_f2(:))   

    % fit_data_gnb = f_get_fit(x_data, y_data, 'gaus_nb');
    % y_fit3 = fit_data_gnb.fit_eq(x_fit, fit_data_gnb.fit_pars);
    % 
    % f1 = fit((1:numel(tun_mean))',tun_mean, 'gauss1');
    % y_fit = f1(x_fit);
    % resid = f1(1:10) - tun_mean;
    % 
    % figure(); hold on;
    % plot(tun_mean, 'k');
    % plot(x_fit, y_fit, '--r');
    % plot(x_fit, y_fit3, '-.g');
    % 
    % 
    % t1 = tinv(0.975,fit_data_gnb.df);
    % 
    % f1.a1 - sqrt(fit_data_gnb.param_cov(1,1))*t1
    % 
    % f1.b1 - sqrt(fit_data_gnb.param_cov(2,2))*t1
    % 
    % f1.c1 - sqrt(fit_data_gnb.param_cov(3,3))*t1


end
% if num_gr > 1
%     if num_reg > 1
%         % fisher lsd t = (mean1 - mean2) / sqrt(MSE(1/n1 - 1/n2))
%         p_all = cell(num_tn,1);
%         tbl_all = cell(num_tn,1);
%         stats_all = cell(num_tn,1);
%         for n_tn = 1:num_tn
%             [p_all{n_tn}, tbl_all{n_tn}, stats_all{n_tn}]  = anova1(cat(1, tun_data_all{:,n_tn}), cat(1, reg_lab_all{:}), 'off');
% 
%             title_tag4 = sprintf('%s; resp counts %s', title_tag3, categories{n_tn});
%             f_dv_plot_anova1(p_all{n_tn}, tbl_all{n_tn}, stats_all{n_tn}, title_tag4);
%         end
%     else
%         [p_all, tbl_all, stats_all]  = anova1(cat(2, tun_data_all{:}),[], 'off');
%         title_tag4 = sprintf('%s; resp counts', title_tag3);
%         f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
%     end
% end

figure; hold on;
for n_reg = 1:num_reg
    has_cells = num_cells_gr2(:,n_reg)>min_cells;
    num_cells4 = num_cells_gr2(has_cells,n_reg);
    tun_data2 = tuned_frac_all2(has_cells,:,:);
    for n_tn = 1:num_tn
        plot(num_cells4, tun_data2(:, n_reg, n_tn), '.', 'color', app.ops.context_types_all_colors2{tn_all(1,n_tn)})
    end
end
title(sprintf('%s; scatter', title_tag3), 'Interpreter', 'none');

%% loco cells

loco_mean = zeros(1, num_reg);
loco_sem = zeros(1, num_reg);
loco_count = zeros(1, num_reg);
loco_data_all = cell(num_reg, 1);
reg_lab_all = cell(num_reg, 1);
for n_reg = 1:num_reg
    has_cells = num_cells_gr2(:,n_reg)>min_cells;
    num_gr_w_cells = sum(has_cells);
    loco_data2 = loco_frac_all2(has_cells,n_reg);
    loco_mean(:,n_reg) = mean(loco_data2,1);
    loco_count(:,n_reg) = sum(has_cells);
    loco_sem(:,n_reg) = std(loco_data2,1)/sqrt(loco_count(:,n_reg) - 1);
    loco_data_all{n_reg} = loco_data2;
    reg_lab_all{n_reg} = repmat(n_reg, num_gr_w_cells, 1);
end

figure; hold on;
bar(categorical(leg_list2,leg_list2), loco_mean, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
for n_br = 1:num_reg
    br1 = bar(n_br, loco_mean(n_br));
    br1.FaceColor = colors1{n_br};
end
if num_gr > 1
    errorbar(1:num_reg, loco_mean, loco_sem, '.k','LineWidth',1);
end
title(sprintf('loco counts; %s', title_tag3), 'Interpreter', 'none');
ylabel('Loco cell fraction');

if num_gr > 1
    [p_loco, tbl_loco, stats_loco]  = anova1(cat(1, loco_data_all{:}), cat(1, reg_lab_all{:}), 'off');
    title_tag4 = sprintf('loco counts; %s', title_tag3);
    f_dv_plot_anova1(p_loco, tbl_loco, stats_loco, title_tag4);
end

%% all imaged neurons

num_cells_all3 = cat(1, num_cells_all{1,:});
num_cells_all4 = sum(num_cells_all3,1);

%if ~app.poolregionsCheckBox.Value
figure; hold on
bar(categorical(leg_list2,leg_list2), num_cells_all4);
title(sprintf('%s %s; cell counts; %d cells total', title_tag, reg_tag, sum(num_cells_all4)), 'Interpreter', 'none');
ylabel('number of cells');
for n_br = 1:numel(num_cells_all4)
    br1 = bar(n_br, num_cells_all4(n_br));
    br1.FaceColor = colors1{n_br};
end
%end

%% all resp cells marg across tn

tun_mean = zeros(1, num_reg);
tun_sem = zeros(1, num_reg);
tun_count = zeros(1, num_reg);
tun_data_all = cell(num_reg, 1);
reg_lab_all = cell(num_reg, 1);
for n_reg = 1:num_reg
    has_cells = num_cells_gr2(:,n_reg)>min_cells;
    num_gr_w_cells = sum(has_cells);
    tun_data2 = tuned_frac_margtn_all2(has_cells,n_reg);
    tun_mean(:,n_reg) = mean(tun_data2,1);
    tun_count(:,n_reg) = sum(has_cells);
    tun_sem(:,n_reg) = std(tun_data2,1)/sqrt(tun_count(:,n_reg) - 1);
    tun_data_all{n_reg} = tun_data2;
    reg_lab_all{n_reg} = repmat(n_reg, num_gr_w_cells, 1);
end

figure; hold on;
bar(categorical(leg_list2,leg_list2), tun_mean, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
for n_br = 1:num_reg
    br1 = bar(n_br, tun_mean(n_br));
    br1.FaceColor = colors1{n_br};
end
if num_gr > 1
    errorbar(1:num_reg, tun_mean, tun_sem, '.k','LineWidth',1);
end
title(sprintf('resp cells marg; %s', title_tag3), 'Interpreter', 'none');
ylabel('Resp cells fraction');

if num_gr > 1
    [p_tun, tbl_tun, stats_tun]  = anova1(cat(1, tun_data_all{:}), cat(1, reg_lab_all{:}), 'off');
    title_tag4 = sprintf('resp cells marg; %s', title_tag3);
    f_dv_plot_anova1(p_tun, tbl_tun, stats_tun, title_tag4, leg_list2);
end

end
