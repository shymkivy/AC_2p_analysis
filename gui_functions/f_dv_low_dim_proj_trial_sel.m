function f_dv_low_dim_proj_trial_sel(app)

method = app.DimredmethodDropDown.Value;
dist_metric = app.LDdistmethodDropDown.Value; % pca isomap
trial_type_sel = app.trialtypeselDropDown.Value;
%trial_type_sel = 'Context';
num_pad = 0;

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

if strcmpi(trial_type_sel, 'Context')
    MMN_idx = 2;
    tn_all = [19, 20];
elseif strcmpi(trial_type_sel, 'Context flip')
    MMN_idx = 1;
    tn_all = [29, 30];
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
    
    tn_all2 = [tn_all, cont_tn_all];
    
    if ~(sum(cont_tn_all > 10) || sum(cont_tn_all < 1))
        [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all2, [], resp_cell_sel);


        cont_tr_all{n_dset} = cont_tn;
        data_all{n_dset} = cat(2,resp_vals{:});
    end
end

data_all2 = cat(1,data_all{:});
cont_tr_all2 = cat(1,cont_tr_all{:});

cont_tn = round(mean(cont_tr_all2));
cont_tn_all = (cont_tn - num_pad):(cont_tn + num_pad);
tn_all2 = [tn_all, cont_tn_all];
cont_idx = find(logical(sum(tn_all2' == cont_tn_all,2)));
ctx_idx = [find(tn_all2 == tn_all(1)), find(tn_all2 == cont_tn), find(tn_all2 == tn_all(2))];

if strcmpi(method, 'pca')

    % reconstruct: data_rec = score*coeff' + mu;
    [coeff,~,~,~,explained,~] = pca(data_all2);
    
    residual_var = round(1 - cumsum(explained/100),4);
    
    lr_data2d = coeff(:,1:2);
    lr_data3d = coeff(:,1:3);

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

figure;
plot([0, 1:numel(residual_var)], [100; residual_var*100], 'o-');
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

figure; hold on
plot3(lr_data3d(cont_idx,1), lr_data3d(cont_idx,2), lr_data3d(cont_idx,3), '.-k')
plot3(lr_data3d(ctx_idx,1), lr_data3d(ctx_idx,2), lr_data3d(ctx_idx,3), '.-k')
for n_tn = 1:numel(tn_all2)
    color1 = app.ops.context_types_all_colors2{tn_all2(n_tn)};
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot3(lr_data3d(n_tn, 1), lr_data3d(n_tn, 2), lr_data3d(n_tn, 3), 'o', 'color', color1, 'LineWidth', 2)
end
title(sprintf('low rank proj trials 3d; %s', title_tag2), 'interpreter', 'none');
grid on

end