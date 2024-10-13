function f_dv_plot_cont_dev_red(app)

bar_sm_ylim = [0 .85];
bar_dist_ylim = [0 7];

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = size(data,1);

[region_num, reg_tag, leg_list, reg_col] = f_dv_get_region_sel_val(app);
num_reg = size(region_num,1);

resp_cell_type = app.ResposivecellstypeDropDown.Value;

if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
    resp_cell_sel = 'All';
else
    resp_cell_sel = 'Resp marg';
end

tn_all = [18 19 20; 28 29 30];
[num_gr, num_tn] = size(tn_all);
num_effdsets = num_dsets*num_gr;

[~, ~, ~, features_proc] = f_dv_get_feature(app, app.plotfeatureDropDown.Value, tn_all, resp_cell_sel);
features_proc2 = reshape(features_proc, num_effdsets, num_reg, num_tn);

fet_all = cat(1,features_proc2{:});
max_pt = round(max(fet_all),2);
min_pt = round(min(fet_all),2);

%%
title_tag2 = sprintf('%s; %s; %s; %s', resp_cell_type, resp_cell_sel, title_tag, reg_tag);
f1 = figure; hold on;
for n_reg = 1:num_reg
    features_proc31 = cat(1,features_proc2{:,n_reg,1});
    features_proc32 = cat(1,features_proc2{:,n_reg,2});
    features_proc33 = cat(1,features_proc2{:,n_reg,3});
    plot3(features_proc31, features_proc32, features_proc33, '.', MarkerSize=1, color=reg_col{n_reg})
end
xlabel('cont'); ylabel('red'); zlabel('dev');
f1.Children.XDir = 'reverse';
f1.Children.YDir = 'reverse';
f1.Children.XLim = [min_pt max_pt];
f1.Children.YLim = [min_pt max_pt];
f1.Children.ZLim = [min_pt max_pt];
grid on
title(sprintf('%s', title_tag2), 'interpreter', 'none')
if num_reg>1
    legend(leg_list);
end
%%

mouse_id = repmat(data.mouse_id, num_gr, 1);
dset_id = num2cell(repmat((1:num_dsets)', num_gr, 1));
[gr_id, gr_tag] = f_dv_combine_data(app, mouse_id, dset_id);
gr_all = unique(gr_id);
num_dgr = numel(gr_all);

plor_pairs = [1, 3; 2, 3; 2, 1];
title_ctx = {'cont', 'red', 'dev'};
num_pl = size(plor_pairs,1);

bar_leg = cell(3,1);
for n_pl = 1:num_pl
    pp2 = plor_pairs(n_pl,:);
    bar_leg{n_pl} = sprintf('%s vs %s', title_ctx{pp2(1)}, title_ctx{pp2(2)});
end

cos_si_all = cell(num_reg, num_pl);
cos_si_mean = zeros(num_reg, num_pl);
dist_all = cell(num_reg, num_pl);
dist_mean = zeros(num_reg, num_pl);
lab_reg = cell(num_reg, num_pl);
lab_cond = cell(num_reg, num_pl);
for n_reg = 1:num_reg
    for n_pl = 1:num_pl
        pp2 = plor_pairs(n_pl,:);
        cos_sim1 = nan(num_dgr,1);
        dist2 = nan(num_dgr,1);
        for n_dgr = 1:num_dgr
            dgr_idx = gr_id == gr_all(n_dgr);
            features_proc31 = cat(1,features_proc2{dgr_idx,n_reg,pp2(1)});
            features_proc32 = cat(1,features_proc2{dgr_idx,n_reg,pp2(2)});
            
            num_cells = numel(features_proc31);
            
            if num_cells > 5
                cos_sim1(n_dgr) = 1 - pdist2(features_proc31', features_proc32', 'cosine');
                dist1 = pdist2(features_proc31', features_proc32', 'euclidean');
                dist2(n_dgr) = dist1/sqrt(num_cells)*sqrt(100);
            end
        end
        cos_sim2 = cos_sim1(~isnan(cos_sim1));
        dist3 = dist2(~isnan(dist2));
        cos_si_all{n_reg, n_pl} = cos_sim2;
        dist_all{n_reg, n_pl} = dist3;
        cos_si_mean(n_reg, n_pl) = mean(cos_sim2, 'omitnan');
        dist_mean(n_reg, n_pl) = mean(dist3, 'omitnan');

        lab_reg{n_reg, n_pl} = repmat(leg_list(n_reg),numel(cos_sim2),1);
        lab_cond{n_reg, n_pl} = repmat(bar_leg(n_pl),numel(cos_sim2),1);
    end
end
%%


for n_pl = 1:num_pl
    f2 = figure; hold on;
    
    pp2 = plor_pairs(n_pl,:);
    for n_reg = 1:num_reg
        features_proc31 = cat(1,features_proc2{:,n_reg,pp2(1)});
        features_proc32 = cat(1,features_proc2{:,n_reg,pp2(2)});

        plot(features_proc31, features_proc32, '.', 'MarkerSize', 1, color=reg_col{n_reg});
        xlabel(title_ctx{pp2(1)}); ylabel(title_ctx{pp2(2)});
        f2.Children.XLim = [min_pt max_pt];
        f2.Children.YLim = [min_pt max_pt];
        title(sprintf('%s; %s; %s; cos SI = %.2f', title_tag2, bar_leg{n_pl}, leg_list{n_reg}, cos_si_mean(n_reg, n_pl)), 'interpreter', 'none')
    end
    if num_reg>1
        legend(leg_list);
    end
end

%%

% for n_reg = 1:num_reg
%     figure;
%     subplot(1,2,1);
%     bar(categorical(categorical(bar_leg)), cos_si_mean(n_reg,:))
%     ylim(bar_sm_ylim);
%     ylabel('Cosine similarity');
%     title('Cosine similarity');
%     subplot(1,2,2);
%     bar(categorical(categorical(bar_leg)), dist_mean(n_reg,:));
%     %ylim(bar_dist_ylim);
%     title('Normalized euclidean distance');
%     ylabel('Normalized euclidean distance');
%     sgtitle(sprintf('%s; %s', title_tag2, leg_list{n_reg}), 'interpreter', 'none');
% end


%%
%labl1 = app.ops.context_types_labels(tn_all(1,:));
%labl2 = app.ops.context_types_labels_trim(tn_all(1,:));

%lab_reg3{n_reg, n_tn} = repmat(leg_list(n_reg), [sum(has_data), 1]);
%lab_tn3{n_reg, n_tn} = repmat(app.ops.context_types_labels_trim(tn_all(1,n_tn)), [sum(has_data), 1]);


cos_si_all2 = cat(1,cos_si_all{:});
dist_all2 = cat(1, dist_all{:});
lab_reg2 = cat(1,lab_reg{:});
lab_cond2 = cat(1,lab_cond{:});

f_plot_bar(cos_si_all, bar_leg, app.ops.context_types_all_colors2(tn_all(1,:)), app.ops.cond_colors);
ylabel('Cosine similarity')
title(title_tag2, 'interpreter', 'none')
% if num_reg>1
%     legend(leg_list);
% end

if app.plotstatsCheckBox.Value
    [p_val,tbl,stats,~] = anovan(cos_si_all2, {lab_reg2, lab_cond2},'varnames', {'Regions', 'Conditions'}, 'model', 'linear');
end

f_plot_bar(dist_all, bar_leg, app.ops.context_types_all_colors2(tn_all(1,:)), app.ops.cond_colors);
ylabel('Normalized euclidean distance')
title(title_tag2, 'interpreter', 'none')

if app.plotstatsCheckBox.Value
    [p_val,tbl,stats,~] = anovan(dist_all2, {lab_reg2, lab_cond2},'varnames', {'Regions', 'Conditions'}, 'model', 'linear');
end

if app.plotstatsCheckBox.Value

    %f_dv_plot_anovan(p_val, tbl, stats, cos_si_all2, {lab_reg2, lab_cond2}, title_tag)
    %fprintf('Regions diff Fanonan(%d)=%.2f, p=%.2e\n', tbl{4, 3}, tbl{2, 6}, p_val(1))
    %fprintf('Trial diff Fanonan(%d)=%.2f, p=%.2e\n', tbl{4, 3}, tbl{3, 6}, p_val(2))

    for n_pl = 1:num_pl
        [p_tun, tbl_tun, stats_tun]  = anova1(cat(1, cos_si_all{:,n_pl}), cat(1, lab_reg{:,n_pl}), 'off');
        title_tag4 = sprintf('cos sim; %s; %s', title_tag2, bar_leg{n_pl});
        f_dv_plot_anova1(p_tun, tbl_tun, stats_tun, title_tag4, leg_list);
    end
    for n_pl = 1:num_pl
        [p_tun, tbl_tun, stats_tun]  = anova1(cat(1, dist_all{:,n_pl}), cat(1, lab_reg{:,n_pl}), 'off');
        title_tag4 = sprintf('euc dist; %s; %s', title_tag2, bar_leg{n_pl});
        f_dv_plot_anova1(p_tun, tbl_tun, stats_tun, title_tag4, leg_list);
    end
end

%%

% if num_reg>1
%     legend(leg_list);
% end

% %%
% cos_si_all = zeros(num_pl,1);
% dist_all = zeros(num_pl,1);
% for n_pl = 1:num_pl
%     pp2 = plor_pairs(n_pl,:);
%     data_all3(or(isnan(data_all2(:,pp2(1))), isnan(data_all2(:,pp2(2)))),:) = [];
%     cos_sim1 = 1 - pdist2(data_all3(:,pp2(1))', data_all3(:,pp2(2))', 'cosine');
%     dist1 = pdist2(data_all3(:,pp2(1))', data_all3(:,pp2(2))', 'euclidean');
% 
%     cos_si_all(n_pl) = cos_sim1;
%     dist_all(n_pl) = dist1;
% end
% 
% %%
% 

% 
% for n_reg = 1:num_reg
%     data_plot = cell(num_dgr, num_tn);
%     for n_tn = 1:num_tn
%         for n_dgr = 1:num_dgr
%             dgr_idx = gr_id == gr_all(n_dgr);
%             data_plot{n_dgr, n_tn} = cat(1,features_proc2{dgr_idx,n_reg,n_tn});
%         end
%     end
% end
% 
% 
% %%
% 
% data_all = cell(num_dsets,num_gr);
% for n_dset = 1:num_dsets
%     data1 = data(n_dset,:);
%     stats1 = cat(1,data1.stats{n_pl});
% 
%     for n_gr = 1:num_gr
%         tn1 = tn_all(n_gr,:);
%         [selected_cells, resp_vals, resp_vals_full] = f_dv_get_resp_vals_cells(app, stats1, tn1, [], resp_cell_sel);
% 
%         resp_vals2 = cat(2,resp_vals{:});
% 
%         num_cells = size(resp_vals2,1);
% 
%         if app.ConverttoZCheckBox.Value
%             st_mean_mean = stats1.stat_trials_mean_mean(selected_cells(:,1));
%             st_mean_sem = stats1.stat_trials_mean_sem(selected_cells(:,1));
%         else
%             st_mean_mean = zeros(num_cells,1);
%             st_mean_sem = ones(num_cells,1);
%         end
% 
%         resp_vals3 = (resp_vals2 - st_mean_mean)./st_mean_sem;
% 
%         data_all{n_dset, n_gr} = resp_vals3;
%     end
% end
% 
% data_all2 = cat(1,data_all{:});
% 
% data_all3 = data_all2;
% 
% max_pt = max(data_all2(:));
% min_pt = round(max_pt*.1,2);
% 
% title_tag2 = sprintf('%s; resp %s; %s', resp_cell_type, resp_cell_sel, title_tag);
% 
% f1 = figure;
% plot3(data_all2(:,1), data_all2(:,2), data_all2(:,3), '.k', 'MarkerSize', 1)
% xlabel('cont'); ylabel('red'); zlabel('dev');
% f1.Children.XDir = 'reverse';
% f1.Children.YDir = 'reverse';
% f1.Children.XLim = [-min_pt max_pt+min_pt];
% f1.Children.YLim = [-min_pt max_pt+min_pt];
% f1.Children.ZLim = [-min_pt max_pt+min_pt];
% grid on
% title(sprintf('%s', title_tag2), 'interpreter', 'none')
% 
% 
% 
% plor_pairs = [1, 3; 2, 3; 2, 1];
% title_ctx = {'cont', 'red', 'dev'};
% 
% num_pl = size(plor_pairs,1);
% 
% cos_si_all = zeros(num_pl,1);
% dist_all = zeros(num_pl,1);
% bar_leg = cell(3,1);
% for n_pl = 1:num_pl
%     pp2 = plor_pairs(n_pl,:);
% 
%     data_all3(or(isnan(data_all2(:,pp2(1))), isnan(data_all2(:,pp2(2)))),:) = [];
%     cos_sim1 = 1 - pdist2(data_all3(:,pp2(1))', data_all3(:,pp2(2))', 'cosine');
%     dist1 = pdist2(data_all3(:,pp2(1))', data_all3(:,pp2(2))', 'euclidean');
% 
%     cos_si_all(n_pl) = cos_sim1;
%     dist_all(n_pl) = dist1;
%     bar_leg{n_pl} = sprintf('%s vs %s', title_ctx{pp2(1)}, title_ctx{pp2(2)});
% 
% 
%     f2 = figure;
%     plot(data_all2(:,pp2(1)), data_all2(:,pp2(2)), '.k', 'MarkerSize', 1);
%     xlabel(title_ctx{pp2(1)}); ylabel(title_ctx{pp2(2)});
%     f2.Children.XLim = [-min_pt max_pt+min_pt];
%     f2.Children.YLim = [-min_pt max_pt+min_pt];
%     title(sprintf('%s; %s, cos SI = %.2f', title_tag2, bar_leg{n_pl}, cos_sim1), 'interpreter', 'none')
% end
% 
% figure;
% subplot(1,2,1);
% bar(categorical(categorical(bar_leg)), cos_si_all)
% ylim(bar_sm_ylim);
% ylabel('Cosine similarity');
% title('Cosine similarity');
% subplot(1,2,2);
% bar(categorical(categorical(bar_leg)), dist_all);
% %ylim(bar_dist_ylim);
% title('Euclidean distance');
% ylabel('Euclidean distance');
% sgtitle(sprintf('%s', title_tag2), 'interpreter', 'none');
% 
% f3 = figure;
% plot(data_all2(:,2), data_all2(:,3), '.k', 'MarkerSize', 1);
% xlabel('red'); ylabel('dev');
% f3.Children.XLim = [-min_pt max_pt+min_pt];
% f3.Children.YLim = [-min_pt max_pt+min_pt];
% title(sprintf('%s; red vs dev', title_tag2), 'interpreter', 'none')
% 
% f4 = figure;
% plot(data_all2(:,2), data_all2(:,1), '.k', 'MarkerSize', 1);
% xlabel('red'); ylabel('cont');
% f4.Children.XLim = [-min_pt max_pt+min_pt];
% f4.Children.YLim = [-min_pt max_pt+min_pt];
% title(sprintf('%s; red vs cont', title_tag2), 'interpreter', 'none')

end