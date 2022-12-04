function f_dv_low_dim_proj_trial_ctx(app)

method = app.DimredmethodDropDown.Value;
dist_metric = app.LDdistmethodDropDown.Value; % pca isomap
num_pad = app.LDcontpadEditField.Value;
trial_type_val = app.trialtypeDropDown.Value;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end
if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
    resp_cell_sel = 'All';
else
    resp_cell_sel = 'Resp marg';
end

num_dsets = size(data,1);

%tn_all = f_dv_get_trial_number(app);

title_tag1 = sprintf('%s; %s', title_tag, method);
if strcmpi(method, 'isomap')
    title_tag1 = sprintf('%s; dist %s', title_tag1, dist_metric);
end

if ~strcmpi(trial_type_val, 'Context') % do context
    fprintf('Running context trials; %d control pad', num_pad);
    MMN_idx = 2;
    tn_all = [19, 20];
    title_tag1 = sprintf('%s; Context', title_tag1);
else
    fprintf('Running context flip trials; %d control pad', num_pad);
    MMN_idx = 1;
    tn_all = [29, 30];
    title_tag1 = sprintf('%s; Context flip', title_tag1);
end

data_all = cell(num_dsets,1);
cont_tr_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    stats1 = cat(1,ddata.stats{n_pl});
    MMN_freq = ddata.MMN_freq{1};
    
    app.ops.context_types_all

    cont_tn = find(app.ops.context_types_all == MMN_freq(MMN_idx));
    cont_tn_all = (cont_tn - num_pad):(cont_tn + num_pad);
    
    tn_all_ctx = [tn_all(1) cont_tn tn_all(2)];
    tn_all2 = [tn_all, cont_tn_all];
    
    if ~(sum(cont_tn_all > 10) || sum(cont_tn_all < 1))
        [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all2, [], resp_cell_sel);


        cont_tr_all{n_dset} = cont_tn;
        data_all{n_dset} = cat(2,resp_vals{:});
    end
end

data_all2 = cat(1,data_all{:});
hasnan1 = logical(sum(isnan(data_all2),2));
data_all2 = data_all2(~hasnan1,:);

cont_tr_all2 = cat(1,cont_tr_all{:});

cont_tn = round(mean(cont_tr_all2));
cont_tn_all = (cont_tn - num_pad):(cont_tn + num_pad);
tn_all2 = [tn_all, cont_tn_all];

cont_idx = find(logical(sum(tn_all2' == cont_tn_all,2)));
ctx_idx = [find(tn_all2 == tn_all(1)), find(tn_all2 == cont_tn), find(tn_all2 == tn_all(2))];

[lr_data2d, lr_data3d, residual_var, residual_var_pca] = f_dv_run_dred(data_all2, method, dist_metric);

title_tag2 = sprintf('%s; resp %s', title_tag1, resp_cell_sel);

figure; hold on;
plot([0, 1:numel(residual_var_pca)], [1; residual_var_pca], 'o-', 'Linewidth', 2);
if ~strcmpi(method, 'pca')
    plot([0, 1:numel(residual_var)], [1; residual_var], 'o-', 'Linewidth', 2);
    legend('PCA', method)
end
title(sprintf('Residual variance proj trials; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

figure; hold on;
plot(lr_data2d(cont_idx,1), lr_data2d(cont_idx,2), '.-k')
plot(lr_data2d(ctx_idx,1), lr_data2d(ctx_idx,2), '.-k')
for n_tn = 1:numel(tn_all2)
    color1 = app.ops.context_types_all_colors2{tn_all2(n_tn)};
    plot(lr_data2d(n_tn, 1), lr_data2d(n_tn, 2), 'o', 'color', color1, 'LineWidth', 2)
end
title(sprintf('low rank proj trials 2d; %s', title_tag2), 'interpreter', 'none');

title_tag3 = sprintf('low rank proj trials; %s', title_tag2);
f_dv_plot3_pc2(lr_data3d, tn_all2, [], title_tag3, app.ops.context_types_all_colors2, 1, {cont_idx, ctx_idx})

end