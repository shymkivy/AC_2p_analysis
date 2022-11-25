function f_dv_low_dim_proj_cell_ctx(app)

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

if ~strcmpi(trial_type_val, 'Context flip') % do context
    fprintf('Running context trials; %d control pad', num_pad);
    MMN_idx = 2;
    tn_all = [19, 20];
else
    fprintf('Running context flip trials; %d control pad', num_pad);
    MMN_idx = 1;
    tn_all = [29, 30];
end

cell_labels = cell(num_dsets,1);
data_all = cell(num_dsets,1);
cont_tr_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    stats1 = cat(1,ddata.stats{n_pl});
    MMN_freq = ddata.MMN_freq{1};
    
    app.ops.context_types_all

    cont_tn = find(app.ops.context_types_all == MMN_freq(MMN_idx));
    cont_tn_all = (cont_tn - num_pad):(cont_tn + num_pad);
    
    tn_all2 = [tn_all, cont_tn_all];
    
    if ~(sum(cont_tn_all > 10) || sum(cont_tn_all < 1))
        [sel_cells, resp_vals, resp_vals2, ~, resp_cells] = f_dv_get_resp_vals_cells(app, stats1, tn_all2, [], resp_cell_sel);
        
        sel_idx = logical(sum(sel_cells,2));
        resp_idx = logical(sum(resp_cells,2));
        [~, max_idx] = max(resp_vals2, [], 2);

        max_idx(~resp_idx) = 0;

        if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
            cell_labels{n_dset} = max_idx;
        else
            cell_labels{n_dset} = max_idx(sel_idx);
        end

        cell_labels{n_dset} = max_idx(sel_idx);
        cont_tr_all{n_dset} = cont_tn;
        data_all{n_dset} = cat(2,resp_vals{:});
    end
end

data_all2 = cat(1,data_all{:})';
cell_labels2 = cat(1,cell_labels{:});
cont_tr_all2 = cat(1,cont_tr_all{:});

cont_tn = round(mean(cont_tr_all2));
cont_tn_all = (cont_tn - num_pad):(cont_tn + num_pad);
tn_all2 = [tn_all, cont_tn_all];
cont_idx = find(logical(sum(tn_all2' == cont_tn_all,2)));
ctx_idx = [find(tn_all2 == tn_all(1)), find(tn_all2 == cont_tn), find(tn_all2 == tn_all(2))];

% reconstruct: data_rec = score*coeff' + mu;
[coeff,~,~,~,explained,~] = pca(data_all2);

residual_var_pca = round(1 - cumsum(explained/100),4);

lr_data2d_pca = coeff(:,1:2);
lr_data3d_pca = coeff(:,1:3);

if strcmpi(method, 'pca')
    lr_data2d = lr_data2d_pca;
    lr_data3d = lr_data3d_pca;

elseif strcmpi(method, 'isomap')
    D = pdist2(data_all2', data_all2', dist_metric); % euclidean, cosine');

    % [Y, R, E] = Isomap(D, n_fcn, n_size, options); 
    %    D = N x N matrix of distances (where N is the number of data points)
    %    n_fcn = neighborhood function ('epsilon' or 'k') 
    %    n_size = neighborhood size (value for epsilon or k)
    options1.display = 0;
    options1.dims = 1:min(numel(tn_all2), 10);
    [Y, R, ~] = IsoMap(D, 'k', 15, options1);
    
    residual_var = R';
    
    %[Y, R, E] = IsoMap(D, 'epsilon', 5);
    
    %[mappedX, mapping] = isomap(data_all2, 2); 
    %figure; plot(mappedX(:,1), mappedX(:,2), 'o')
    
    lr_data2d = Y.coords{2}';
    lr_data3d = Y.coords{3}';
    
    title_tag1 = sprintf('%s; dist %s', title_tag1, dist_metric);
end

title_tag2 = sprintf('%s; resp %s', title_tag1, resp_cell_sel);

figure; hold on;
plot([0, 1:numel(residual_var_pca)], [100; residual_var_pca*100], 'o-', 'LineWidth', 2);
if ~strcmpi(method, 'pca')
    plot([0, 1:numel(residual_var)], [100; residual_var*100], 'o-', 'LineWidth', 2);
    legend('PCA', method);
end
title(sprintf('Residual variance proj trials; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

%%
figure; hold on
plot(lr_data2d(:,1), lr_data2d(:,2), 'ok', 'LineWidth', 1)
for n_cell = 1:numel(cell_labels2)
    if cell_labels2(n_cell)
        color1 = app.ops.context_types_all_colors2{tn_all2(cell_labels2(n_cell))};
    else
        color1 = [.6 .6 .6];
    end
    plot(lr_data2d(n_cell, 1), lr_data2d(n_cell, 2), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
end
title(sprintf('low rank proj cells 2d; %s',title_tag2), 'interpreter', 'none');

figure; hold on
plot3(lr_data3d(:,1), lr_data3d(:,2), lr_data3d(:,3), 'ok', 'LineWidth', 1)
for n_cell = 1:numel(cell_labels2)
    if cell_labels2(n_cell)
        color1 = app.ops.context_types_all_colors2{tn_all2(cell_labels2(n_cell))};
    else
        color1 = [.6 .6 .6];
    end
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot3(lr_data3d(n_cell, 1), lr_data3d(n_cell, 2), lr_data3d(n_cell, 3), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
end
title(sprintf('low rank proj cells 3d; %s', title_tag2), 'interpreter', 'none');
grid on

end