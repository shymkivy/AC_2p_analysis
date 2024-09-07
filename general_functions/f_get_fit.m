function data = f_get_fit(data_x, data_y, type, init_params)
% to get std of the parameters, we need to compute Jacobian matrix, and
% inv(J'*J)*sigma^2 = covariance matrix wrt params
% first col is deriv of f wrt b0, eval at all x and computed param vals
% second col is deric of f wrt b1


err = @(y, x) sum((y - x).^2);
n = numel(data_x);

if strcmpi(type, 'linear')
    if ~exist('init_params', 'var')
        init_params = [mean(data_y)-min(data_y), 0];
    end
    fit_eq = @(x, pars) pars(1,:) + pars(2,:).*x;
    k = 1; % free params
elseif strcmpi(type, 'exp')
    if ~exist('init_params', 'var')
        init_params = [max(data_y), 1, min(data_y)];
    end
    fit_eq = @(x, pars) pars(1,:).*exp(-1./pars(2,:).*x) + pars(3,:);
    k=2;
elseif strcmpi(type, 'gaus_nb')
    if ~exist('init_params', 'var')
        init_params = [max(data_y)-min(data_y), mean(data_x), 1];
    end
    fit_eq = @(x, pars) pars(1,:).*exp(-(((x-pars(2,:))./pars(3,:)).^2)); % pars(4)*x + pars(5)
    k=2;
elseif strcmpi(type, 'gaus')
    if ~exist('init_params', 'var')
        init_params = [max(data_y)-min(data_y), mean(data_x), 1, min(data_y)];
    end
    fit_eq = @(x, pars) pars(1,:).*exp(-(((x-pars(2,:))./pars(3,:)).^2)) +pars(4,:); % pars(4)*x + pars(5)
    k=3;
elseif strcmpi(type, 'gaus+line')
    if ~exist('init_params', 'var')
        init_params = [max(data_y)-min(data_y), mean(data_x), 1, min(data_y), 0];
    end
    fit_eq = @(x, pars) pars(1,:).*exp(-(((x-pars(2,:))./pars(3,:)).^2)) +pars(4,:) + pars(5,:).*x; % pars(4)*x + pars(5)
    k=4;
end

df = n-k-1;

functionToMinimize = @(pars, x, y)(norm(fit_eq(x, pars) - y));
targetFunctionForFminseardch = @(pars)(functionToMinimize(pars, data_x, data_y));
options = optimset('MaxIter', 1000*(k+1));
minPars = fminsearch(targetFunctionForFminseardch, init_params(:), options);

resid = fit_eq(data_x, minPars) - data_y;

err_total = err(data_y, mean(data_y));
err_resid = err(data_y, fit_eq(data_x, minPars));

r_sq = 1 - err_resid/err_total;
r_sq_adj = 1 - (1-r_sq) * (n - 1)/(df); % (n - 1)/(n-p-1) n is number data, p is number params excluding constant term

% calculate parameter covariance
h = 1e-5*ones(size(minPars'));

J = @(p,h,F, x)(F(x, repmat(p,size(p'))+diag(h))-F(x, repmat(p, size(p'))))./h;

J_data = J(minPars, h, fit_eq, data_x);

var_er = sum(resid.^2)/(df);

largesigma = inv(J_data'*J_data)*var_er;

data.df = df;
data.fit_type = type;
data.fit_eq = fit_eq;
data.fit_eq_str = func2str(fit_eq);
data.fit_pars = minPars;
data.residuals = resid;
data.err_total = err_total;
data.err_resid = err_resid;
data.r_sq = r_sq;
data.r_sq_adj = r_sq_adj;
data.param_cov = largesigma;


 % linear regression math way
% data_x_mean2 = [ones(10,1), data_x_mean'];
% 
% data_x_mean2*minPars'
% 
% X = data_x_mean2
% 
% beta = inv(X'*X)*X'*data_y_mean';
% 
% y_pred = data_x_mean2*beta;
% 
% y_res = data_y_mean' - y_pred;
% 
% sigma = sqrt(mean(y_res.^2));
% 
% sigma_cov = sigma.^2*inv(X'*X);

end