function [lr_data, residual_var, residual_var_pca, subtr_mean] = f_dv_run_dred(data_all2, params)

max_comp = min(min(size(data_all2)), 10);

if params.subtract_mean
    subtr_mean = mean(data_all2,1);
    data_all2 = data_all2 - subtr_mean;
else
    subtr_mean = zeros(max_comp,1);
end

% reconstruct: data_rec = score*coeff' + mu;
[coeff,score,~,~,explained,mu] = pca(data_all2);

SS_pca = diag(score'*score);

%datab = score*coeff';

residual_var_pca = round(1 - cumsum(explained/100),4);

if params.scale_by_var
    coeff = coeff*diag(sqrt(SS_pca));
end

lr_data_pca = coeff;

if strcmpi(params.method, 'pca')
    lr_data = lr_data_pca;
    
    residual_var = residual_var_pca;
elseif strcmpi(params.method, 'svd')
    [U,S,V] = svd(data_all2);
    
    %datab = U * S * V';
    if params.scale_by_var
        V = V*diag(diag(S));
    end

    lr_data = V;
    
    SS = diag(S).^2;
    explained = SS/sum(SS);
    residual_var = round(1 - cumsum(explained),4);
    
elseif strcmpi(params.method, 'isomap')
    D = pdist2(data_all2', data_all2', params.dist_metric); % euclidean, cosine');

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
    
    lr_data = Y.coords{end}';
end

end