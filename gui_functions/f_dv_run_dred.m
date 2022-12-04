function [lr_data2d, lr_data3d, residual_var, residual_var_pca] = f_dv_run_dred(data_all2, method, dist_metric)

max_comp = min(min(size(data_all2)), 10);

% reconstruct: data_rec = score*coeff' + mu;
[coeff,~,~,~,explained,~] = pca(data_all2);

residual_var_pca = round(1 - cumsum(explained/100),4);

lr_data2d_pca = coeff(:,1:2);
lr_data3d_pca = coeff(:,1:3);

if strcmpi(method, 'pca')
    lr_data2d = lr_data2d_pca;
    lr_data3d = lr_data3d_pca;
    
    residual_var = residual_var_pca;
    
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
end

end