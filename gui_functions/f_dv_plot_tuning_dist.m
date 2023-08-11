function f_dv_plot_tuning_dist(app)

plot_ens = app.plotensamblesCheckBox.Value;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

tn_all = f_dv_get_trial_number(app);
[num_gr, num_tn] = size(tn_all);
num_dsets = numel(data.experiment);
reg_all = app.ops.regions_to_analyze;

[region_num, reg_tag] = f_dv_get_region_sel_val(app);
num_regions = numel(region_num);

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
    
    loco_cells = cat(1,[stats1.loco_resp_cells])';
    
        %reg_cell_idx = logical(sum(reg_cell_labels == region_num,2));
    for n_gr = 1:num_gr
        tn1 = tn_all(n_gr, :);
        mouse_id{n_gr, n_dset} = data1.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;
        
        [~, ~, ~, ~, resp_cells] = f_dv_get_resp_vals_cells(app, stats1, tn1);
        
        cell_is_resp_reg = zeros(num_cells, num_regions, num_tn);
        loco_cells_reg = zeros(num_cells, num_regions);
        num_cells1 = zeros(1, num_regions);
        for n_reg = 1:num_regions
            reg_idx = reg_cell_labels == region_num(n_reg);
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
tuned_frac_marg_all = cell(num_gr, 1);
tuned_frac_margtn_all = cell(num_gr, 1);
loco_frac_all = cell(num_gr, 1);
loco_frac_marg_all = cell(num_gr, 1);
tuned_frac_margtn_marg_all = cell(num_gr, 1);
num_cells_gr = cell(num_gr, 1);

for n_gr = 1:num_gr
    idx1 = gr_id_uq(n_gr) == gr_id;
    resp_cell_all3 = cat(1, resp_cell_all2{idx1});
    num_cells_all3 = sum(cat(1, num_cells_all2{idx1}),1);
    loco_cells_all3 = cat(1, loco_cells_all2{idx1});
    
    num_cells_gr{n_gr} = num_cells_all3;
    tuned_frac_all{n_gr} = sum(resp_cell_all3,1)./num_cells_all3;
    tuned_frac_marg_all{n_gr} = sum(logical(sum(resp_cell_all3,2)),1)/sum(num_cells_all3);
    
    tuned_frac_margtn_all{n_gr} = sum(logical(sum(resp_cell_all3,3)),1)/sum(num_cells_all3);
    tuned_frac_margtn_marg_all{n_gr} = sum(logical(sum(logical(sum(resp_cell_all3,3)),2)),1)/sum(num_cells_all3);
    
    loco_frac_all{n_gr} = sum(loco_cells_all3,1)./num_cells_all3;
    loco_frac_marg_all{n_gr} = sum(logical(sum(loco_cells_all3,2)),1)/sum(num_cells_all3);
end

title_tag3 = sprintf('%s; stats by %s, %d groups', title_tag2, gr_tag, num_gr);

num_cells_gr2 = cat(1, num_cells_gr{:});
tuned_frac_all2 = cat(1, tuned_frac_all{:});
tuned_frac_marg_all2 = cat(1, tuned_frac_marg_all{:});
loco_frac_all2 = cat(1, loco_frac_all{:});
loco_frac_marg_all2 = cat(1, loco_frac_marg_all{:});

tuned_frac_margtn_all2 = cat(1, tuned_frac_margtn_all{:});
tuned_frac_margtn_marg_all2 = cat(1, tuned_frac_margtn_marg_all{:});

if app.poolregionsCheckBox.Value
    num_regions2 = 1;
    num_cells_gr3 = sum(num_cells_gr2,2);
    tun_data = tuned_frac_marg_all2;
    loco_data = loco_frac_marg_all2;
    tun_marg_data = tuned_frac_margtn_marg_all2;
    reg_all2 = {'Comb reg'};
    colors1 = {[.6 .6 .6]};
else
    num_regions2 = num_regions;
    num_cells_gr3 = num_cells_gr2;
    tun_data = tuned_frac_all2;
    loco_data = loco_frac_all2;
    tun_marg_data = tuned_frac_margtn_all2;
    reg_all2 = reg_all;
    colors1 = app.ops.cond_colors;
end

min_cells = 10;

tun_mean = zeros(num_tn,num_regions2);
tun_sem = zeros(num_tn,num_regions2);
tun_count = zeros(1, num_regions2);
tun_data_all = cell(num_regions2,num_tn);
reg_lab_all = cell(num_regions2, 1);
for n_reg = 1:num_regions2
    has_cells = num_cells_gr3(:,n_reg)>min_cells;
    num_cells_tn = sum(has_cells);
    tun_data2 = tun_data(has_cells,n_reg,:);
    tun_mean(:,n_reg) = mean(tun_data2,1);
    tun_count(:,n_reg) = sum(has_cells);
    tun_sem(:,n_reg) = std(tun_data2,1)/sqrt(tun_count(:,n_reg) - 1);
    for n_tn = 1:num_tn
        tun_data_all{n_reg, n_tn} = tun_data2(:,1,n_tn);
    end
    reg_lab_all{n_reg} = repmat(n_reg, num_cells_tn, 1);
end

figure; hold on;
b1 = bar(categorical(categories,categories), tun_mean, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
for n_bar = 1:num_regions2
    b1(n_bar).FaceColor = colors1{n_bar};
    b1(n_bar).EdgeColor = 'k';
    b1(n_bar).LineWidth = 1;
end
xb = nan(num_tn, num_regions2);
if num_gr > 1
    for n_bar = 1:num_regions2
        xb(:, n_bar) = b1(n_bar).XEndPoints;
    end
    errorbar(xb, tun_mean, tun_sem, '.k','LineWidth',1);
end
title(sprintf('%s; resp counts', title_tag3), 'Interpreter', 'none');
ylabel('Resp cell fraction');

if num_gr > 1
    if num_regions2 > 1
        % fisher lsd t = (mean1 - mean2) / sqrt(MSE(1/n1 - 1/n2))
        p_all = cell(num_tn,1);
        tbl_all = cell(num_tn,1);
        stats_all = cell(num_tn,1);
        for n_tn = 1:num_tn
            [p_all{n_tn}, tbl_all{n_tn}, stats_all{n_tn}]  = anova1(cat(1, tun_data_all{:,n_tn}), cat(1, reg_lab_all{:}), 'off');

            title_tag4 = sprintf('%s; resp counts %s', title_tag3, categories{n_tn});
            f_dv_plot_anova1(p_all{n_tn}, tbl_all{n_tn}, stats_all{n_tn}, title_tag4);
        end
    else
        [p_all, tbl_all, stats_all]  = anova1(cat(2, tun_data_all{:}),[], 'off');
        title_tag4 = sprintf('%s; resp counts', title_tag3);
        f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
    end
end

figure; hold on;
for n_reg = 1:num_regions2
    has_cells = num_cells_gr3(:,n_reg)>min_cells;
    num_cells4 = num_cells_gr3(has_cells,n_reg);
    tun_data2 = tun_data(has_cells,:,:);
    for n_tn = 1:num_tn
        plot(num_cells4, tun_data2(:, n_reg, n_tn), '.', 'color', app.ops.context_types_all_colors2{tn_all(1,n_tn)})
    end
end
title(sprintf('%s; scatter', title_tag3), 'Interpreter', 'none');

%% loco cells

loco_mean = zeros(1, num_regions2);
loco_sem = zeros(1, num_regions2);
loco_count = zeros(1, num_regions2);
loco_data_all = cell(num_regions2, 1);
reg_lab_all = cell(num_regions2, 1);
for n_reg = 1:num_regions2
    has_cells = num_cells_gr3(:,n_reg)>min_cells;
    num_cells_tn = sum(has_cells);
    loco_data2 = loco_data(has_cells,n_reg);
    loco_mean(:,n_reg) = mean(loco_data2,1);
    loco_count(:,n_reg) = sum(has_cells);
    loco_sem(:,n_reg) = std(loco_data2,1)/sqrt(loco_count(:,n_reg) - 1);
    loco_data_all{n_reg} = loco_data2;
    reg_lab_all{n_reg} = repmat(n_reg, num_cells_tn, 1);
end

figure; hold on;
bar(categorical(reg_all2,reg_all2), loco_mean, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
for n_br = 1:num_regions2
    br1 = bar(n_br, loco_mean(n_br));
    br1.FaceColor = colors1{n_br};
end
if num_gr > 1
    errorbar(1:num_regions2, loco_mean, loco_sem, '.k','LineWidth',1);
end
title(sprintf('%s; loco counts', title_tag3), 'Interpreter', 'none');
ylabel('Loco cell fraction');

if num_gr > 1
    [p_loco, tbl_loco, stats_loco]  = anova1(cat(1, loco_data_all{:}), cat(1, reg_lab_all{:}), 'off');
    title_tag4 = sprintf('%s; loco counts', title_tag3);
    f_dv_plot_anova1(p_loco, tbl_loco, stats_loco, title_tag4);
end

%% all imaged neurons

num_cells_all3 = cat(1, num_cells_all{1,:});
num_cells_all4 = sum(num_cells_all3,1);

%if ~app.poolregionsCheckBox.Value
figure; hold on
bar(categorical(reg_all,reg_all), num_cells_all4);
title(sprintf('%s %s; cell counts; %d cells total', title_tag, reg_tag, sum(num_cells_all4)), 'Interpreter', 'none');
ylabel('number of cells');
for n_br = 1:numel(num_cells_all4)
    br1 = bar(n_br, num_cells_all4(n_br));
    br1.FaceColor = app.ops.cond_colors{n_br};
end
%end

%% all resp cells marg across tn

tun_mean = zeros(1, num_regions2);
tun_sem = zeros(1, num_regions2);
tun_count = zeros(1, num_regions2);
tun_data_all = cell(num_regions2, 1);
reg_lab_all = cell(num_regions2, 1);
for n_reg = 1:num_regions2
    has_cells = num_cells_gr3(:,n_reg)>min_cells;
    num_cells_tn = sum(has_cells);
    tun_data2 = tun_marg_data(has_cells,n_reg);
    tun_mean(:,n_reg) = mean(tun_data2,1);
    tun_count(:,n_reg) = sum(has_cells);
    tun_sem(:,n_reg) = std(tun_data2,1)/sqrt(tun_count(:,n_reg) - 1);
    tun_data_all{n_reg} = tun_data2;
    reg_lab_all{n_reg} = repmat(n_reg, num_cells_tn, 1);
end

figure; hold on;
bar(categorical(reg_all2,reg_all2), tun_mean, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
for n_br = 1:num_regions2
    br1 = bar(n_br, tun_mean(n_br));
    br1.FaceColor = colors1{n_br};
end
if num_gr > 1
    errorbar(1:num_regions2, tun_mean, tun_sem, '.k','LineWidth',1);
end
title(sprintf('%s; resp cells marg', title_tag3), 'Interpreter', 'none');
ylabel('Resp cells fraction');


if num_gr > 1
    [p_tun, tbl_tun, stats_tun]  = anova1(cat(1, tun_data_all{:}), cat(1, reg_lab_all{:}), 'off');
    title_tag4 = sprintf('%s; resp cells marg', title_tag3);
    f_dv_plot_anova1(p_tun, tbl_tun, stats_tun, title_tag4);
end

end
