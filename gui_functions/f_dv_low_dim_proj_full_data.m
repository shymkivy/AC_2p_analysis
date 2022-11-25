function f_dv_low_dim_proj_full_data(app)

method = app.DimredmethodDropDown.Value;
dist_metric = app.LDdistmethodDropDown.Value;

cdata = f_dv_get_cdata(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([app.ddata.num_planes]);
end

stats1 = cat(1,app.ddata.stats{n_pl});

tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';
stim_times = app.ddata.stim_frame_index{1};
mmn_freq = app.ddata.MMN_freq{1};
trial_types = app.ddata.trial_types{1};

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[~, trial_frames] = f_dv_compute_window_t(trial_window, cdata.volume_period);

firing_rate = cat(1,cdata.S_sm);
firing_rate = f_normalize(firing_rate, 'norm_mean_std');

[sel_cells, ~, resp_vals2, ~, resp_cells] = f_dv_get_resp_vals_cells(app, stats1, tn_all);

sel_idx = logical(sum(sel_cells,2));
resp_idx = logical(sum(resp_cells,2));
[~, resp_type] = max(resp_vals2, [], 2);
resp_type(~resp_idx) = 0;
resp_type2 = resp_type(sel_idx);

firing_rate2 = firing_rate(sel_idx,:);
num_cells = size(firing_rate2,1);

trial_data_sort = f_get_stim_trig_resp(firing_rate2, stim_times, trial_frames);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

trial_idx = logical(sum(tt_all == trial_types_wctx,2));
trial_data_sort2 = trial_data_sort_wctx(:,:,trial_idx);

firing_rate3 = reshape(trial_data_sort2, num_cells, []);

% remove inactive cells
active_cells = sum(firing_rate3,2) ~= 0;
firing_rate3(~active_cells,:) = [];
resp_type2(~active_cells) = [];

title_tag1 = sprintf('%s', method);

if strcmpi(method, 'pca')

    % reconstruct: data_rec = score*coeff' + mu;
    [coeff,~,~,~,explained,~] = pca(firing_rate3');
    
    residual_var = round(1 - cumsum(explained/100),4);
    
    lr_data2d = coeff(:,1:2);
    lr_data3d = coeff(:,1:3);

elseif strcmpi(method, 'isomap')
    D = pdist2(firing_rate3, firing_rate3, dist_metric);

    % [Y, R, E] = Isomap(D, n_fcn, n_size, options); 
    %    D = N x N matrix of distances (where N is the number of data points)
    %    n_fcn = neighborhood function ('epsilon' or 'k') 
    %    n_size = neighborhood size (value for epsilon or k)
    options1.display = 0;
    [Y, R, ~] = IsoMap(D, 'k', 15, options1);
    
    residual_var = R';
    
    %[Y, R, E] = IsoMap(D, 'epsilon', 1);
    
    %[mappedX, mapping] = isomap(data_all2, 2); 
    %figure; plot(mappedX(:,1), mappedX(:,2), 'o')
    
    lr_data2d = Y.coords{2}';
    lr_data3d = Y.coords{3}';
    
    title_tag1 = sprintf('%s; dist %s', title_tag1, dist_metric);
end

title_tag2 = sprintf('%s; resp %s; %d cells', title_tag1, app.ResponsivecellsselectDropDown.Value, num_cells);

figure;
plot([0, 1:10], [1; residual_var], 'o-');
title(sprintf('Residual variance proj cells; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

figure; hold on
plot(lr_data2d(:,1), lr_data2d(:,2), 'ok', 'LineWidth', 1)
for n_cell = 1:num_cells
    if resp_type2(n_cell)
        color1 = app.ops.context_types_all_colors2{tn_all(resp_type2(n_cell))};
    else
        color1 = [.6 .6 .6];
    end
    plot(lr_data2d(n_cell, 1), lr_data2d(n_cell, 2), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
end
title(sprintf('low rank proj cells 2d; %s',title_tag2), 'interpreter', 'none');

figure; hold on
plot3(lr_data3d(:,1), lr_data3d(:,2), lr_data3d(:,3), 'ok', 'LineWidth', 1)
for n_cell = 1:num_cells
    if resp_type2(n_cell)
        color1 = app.ops.context_types_all_colors2{tn_all(resp_type2(n_cell))};
    else
        color1 = [.6 .6 .6];
    end
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot3(lr_data3d(n_cell, 1), lr_data3d(n_cell, 2), lr_data3d(n_cell, 3), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
end
title(sprintf('low rank proj cells 3d; %s', title_tag2), 'interpreter', 'none');
grid on



% 
% no_dims = round(intrinsic_dim(firing_rate2, 'GMST'));
% 
% [X, labels] = generate_data('helix', 2000);
% figure, scatter3(X(:,1), X(:,2), X(:,3), 5, labels); title('Original dataset'), drawnow
% no_dims = round(intrinsic_dim(X, 'MLE'));
% disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]);
% 
% no_dims = 2
% 
% [mappedX, mapping] = compute_mapping(firing_rate2, 'PCA', 3);	
% figure, scatter3(mappedX(:,1), mappedX(:,2), mappedX(:,3), 5); title('Result of PCA');
% 
% [mappedX, mapping] = compute_mapping(firing_rate2', 'Isomap', 2);	
% figure, scatter(mappedX(:,1), mappedX(:,2), 5); title('Result of Isomap');
% 
% [mappedX, mapping] = compute_mapping(firing_rate2, 'GPLVM', no_dims, 7);	
% figure, scatter(mappedX(:,1), mappedX(:,2), 5, labels); title('Result of GPLVM'); drawnow % (mapping.conn_comp)
% 

end