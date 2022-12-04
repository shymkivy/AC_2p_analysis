function [lr_data, residual_var, residual_var_pca] = f_dv_run_dred2(data_all2, method, dist_metric)

max_comp = min(min(size(data_all2)), 10);

% reconstruct: data_rec = score*coeff' + mu;
[coeff,~,~,~,explained,~] = pca(data_all2);

max_comp = min(min(numel(explained)), 10);


residual_var_pca = round(1 - cumsum(explained/100),4);

lr_data = cell(max_comp,1);
if strcmpi(method, 'pca')
    for n_comp = 1:max_comp
        lr_data{n_comp} = coeff(:,1:n_comp);
    end
    
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
    
    for n_comp = 1:max_comp
        lr_data{n_comp} = Y.coords{n_comp}';
    end
end

end