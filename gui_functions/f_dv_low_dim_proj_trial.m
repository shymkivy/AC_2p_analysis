function f_dv_low_dim_proj_trial(app)

method = app.DimredmethodDropDown.Value;
dist_metric = app.LDdistmethodDropDown.Value; % pca isomap

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = size(data,1);

tn_all = f_dv_get_trial_number(app);

title_tag1 = sprintf('%s; %s', title_tag, method);

data_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    stats1 = cat(1,data(n_dset,:).stats{n_pl});

    if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
        resp_cell_sel = 'All';
    else
        resp_cell_sel = 'Resp marg';
    end

    [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], resp_cell_sel);

    data_all{n_dset} = cat(2,resp_vals{:});
end

data_all2 = cat(1,data_all{:});
hasnan1 = logical(sum(isnan(data_all2),2));
data_all2 = data_all2(~hasnan1,:);
max_comp = min(numel(tn_all), 10);

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
    options1.dims = 1:max_comp;
    [Y, R, ~] = IsoMap(D, 'k', 15, options1);
    
    residual_var = R';
    
    %[Y, R, E] = IsoMap(D, 'epsilon', 1);
    
    %[mappedX, mapping] = isomap(data_all2, 2); 
    %figure; plot(mappedX(:,1), mappedX(:,2), 'o')
    
    lr_data2d = Y.coords{2}';
    lr_data3d = Y.coords{3}';
    
    title_tag1 = sprintf('%s; dist %s', title_tag1, dist_metric);
end

title_tag2 = sprintf('%s; resp %s', title_tag1, resp_cell_sel);

figure; hold on;
plot([0, 1:numel(residual_var_pca)], [1; residual_var_pca], 'o-', 'Linewidth', 2);
if ~strcmpi(method, 'pca')
    plot([0, 1:numel(residual_var_pca)], [1; residual_var], 'o-', 'Linewidth', 2);
    legend('PCA', method)
end
title(sprintf('Residual variance proj trials; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

figure; hold on
plot(lr_data2d(:,1), lr_data2d(:,2), '.-k')
for n_tn = 1:numel(tn_all)
    plot(lr_data2d(n_tn, 1), lr_data2d(n_tn, 2), 'o', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
end
title(sprintf('low rank proj trials 2d; %s', title_tag2), 'interpreter', 'none');

figure; hold on
plot3(lr_data3d(:,1), lr_data3d(:,2), lr_data3d(:,3), '.-k')
for n_tn = 1:numel(tn_all)
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot3(lr_data3d(n_tn, 1), lr_data3d(n_tn, 2), lr_data3d(n_tn, 3), 'o', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
end
title(sprintf('low rank proj trials 3d; %s', title_tag2), 'interpreter', 'none');
grid on

end