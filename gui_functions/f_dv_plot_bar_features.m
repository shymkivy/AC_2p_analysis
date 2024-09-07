function f_dv_plot_bar_features(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

tn_all = f_dv_get_trial_number(app);
[num_gr, num_tn] = size(tn_all);

transp = app.StimtranspEditField.Value;
freq_col = app.stimcolorSpinner.Value;

[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val2(app);
num_reg = size(region_num,1);
num_dsets = size(data,1);

title_tag3 = sprintf('%s; %s; reg %s; %s', title_tag, app.ResponsivecellsselectDropDown.Value, reg_tag, app.plotfeatureDropDown.Value);

features1 = cell(num_dsets, num_gr, num_reg, num_tn);
num_cells_all = zeros(num_dsets, num_gr, num_reg);

for n_gr = 1:num_gr
    [features0, sel_cells0, area_labels0] = f_dv_get_feature2(app, app.plotfeatureDropDown.Value, data, tn_all(n_gr,:), n_pl);
    for n_dset = 1:num_dsets
        for n_reg = 1:num_reg
            reg_idx1 = logical(sum(area_labels0{n_dset} == region_num(n_reg,:),2));
            
            if sum(reg_idx1)
                for n_tn = 1:num_tn
                    feat1 = features0{n_dset, n_tn};
                    sel1 = sel_cells0{n_dset, n_tn};
                    
                    feat2 = feat1(reg_idx1);
                    sel2 = sel1(reg_idx1);

                    feat3 = feat2(sel2);
                    
                    features1{n_dset, n_gr, n_reg, n_tn} = feat3;
                    num_cells_all(n_dset, n_gr, n_reg, n_tn) = sum(sel2);
                end
            end
        end
    end
end

num_effdsets = num_dsets*num_gr;
mouse_id = cell(num_gr, num_dsets);
dset_id = zeros(num_gr, num_dsets);
for n_dset = 1:num_dsets
    data1 =  data(n_dset,:);
    for n_gr = 1:num_gr
        mouse_id{n_gr, n_dset} = data1.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;
    end
end
mouse_id2 = reshape(mouse_id, num_effdsets, 1);
dset_id2 = reshape(dset_id, num_effdsets, 1);
[sd_id, sd_tag] = f_dv_combine_data(app, mouse_id2, dset_id2);
sd_all = unique(sd_id);
num_sd = numel(sd_all);

features2 = reshape(features1, num_effdsets, num_reg, num_tn);

features3 = cell(num_reg, num_tn);
num_cells_all3 = cell(num_reg, num_tn);
lab_reg3 = cell(num_reg, num_tn);
lab_tn3 = cell(num_reg, num_tn);
for n_reg = 1:num_reg
    for n_tn = 1:num_tn
        if num_sd == 1
            feat_32 = cat(1,features2{:,n_reg,n_tn});
            num_cells_32 = ones(numel(feat_32),1);
            has_data = true(numel(feat_32),1);
        else
            feat_32 = zeros(num_sd,1);
            num_cells_32 = zeros(num_sd,1);
            has_data = false(num_sd,1);
            for n_sd = 1:num_sd
                idx1 = sd_all(n_sd) == sd_id;
                
                feat2 = cat(1,features2{idx1,n_reg,n_tn});
                
                num_cells333 = numel(feat2);
                if num_cells333
                    feat_32(n_sd) = mean(feat2);
                    num_cells_32(n_sd) = num_cells333;
                    has_data(n_sd) = 1;
                end
            end
        end
        features3{n_reg,n_tn} = feat_32(has_data);
        num_cells_all3{n_reg,n_tn} = num_cells_32(has_data);
        
        
        lab_reg3{n_reg,n_tn} = repmat(leg_list(n_reg), [sum(has_data), 1]);
        lab_tn3{n_reg,n_tn} = repmat(app.ops.context_types_labels_trim(tn_all(1,n_tn)), [sum(has_data), 1]);
    end
end

labl1 = app.ops.context_types_labels(tn_all(1,:));
labl2 = app.ops.context_types_labels_trim(tn_all(1,:));

f_plot_bar(features3, labl2, app.ops.context_types_all_colors2(tn_all(1,:)), app.ops.cond_colors);
ylabel(app.plotfeatureDropDown.Value)
title(title_tag3, 'interpreter', 'none')

feat4 = cat(1,features3{:});
lab_reg4 = cat(1,lab_reg3{:});
lab_tn4 = cat(1,lab_tn3{:});

[p_val,tbl,stats,~] = anovan(feat4, {lab_reg4, lab_tn4},'varnames', {'Regions', 'Trials'}, 'model', 'linear');

f_dv_plot_anovan(p_val, tbl, stats, feat4, {lab_reg4, lab_tn4}, title_tag)
fprintf('Regions diff Fanonan(%d)=%.2f, p=%.2e\n', tbl{4, 3}, tbl{2, 6}, p_val(1))
fprintf('Trial diff Fanonan(%d)=%.2f, p=%.2e\n', tbl{4, 3}, tbl{3, 6}, p_val(2))

for n_tn = 1:num_tn
    [p_all,tbl_all,stats_all] = anova1(cat(1,features3{:,n_tn}), cat(1,lab_reg3{:,n_tn}), 'off');
    title_tag5 = sprintf('between reg %s; %s; %s', labl1{n_tn}, title_tag3); 
    f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag5, leg_list);
end


% if ~app.MarginalizedistCheckBox.Value
%     [p_all, tbl_all, stats_all]  = anova1(cat(1, feat_pool{:}),cat(1, lab_pool{:}), 'off');
%     title_tag4 = sprintf('%s; stats', title_tag2);
%     f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
% end

end