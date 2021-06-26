function data_dim_est = f_ensemble_comp_data_dim2(firing_rate, params)
% input either 3D tials data (Cell x Time X Trial)
%           or 2D trial data (Cell x Trial)
%% parameters
if ~exist('params', 'var') || isempty(params)
    params = struct;
end

normalize1 = f_get_param(params, 'normalize', 'norm_mean_std');  % 'norm_mean_std', 'norm_mean' 'none'
shuffle_method = f_get_param(params, 'shuffle_method', 'circ_shift');     % 'circ_shift' or 'scramble'
total_dim_thresh = f_get_param(params, 'total_dim_thresh', .7);
dim_est_num_reps = f_get_param(params, 'dim_est_num_reps', 50);
%corr_comp_thresh = f_get_param(params, 'corr_comp_thresh', .90);
plot_stuff = f_get_param(params, 'plot_stuff', 0);

%%
ndims1 = ndims(firing_rate);
if ndims1 == 3
    num_cells = size(firing_rate,1);
    firing_rate = reshape(firing_rate, num_cells,[]);
end  

% active_cells = sum(firing_rate,2) ~= 0;
% firing_rate(~active_cells,:) = [];
% ens_out.active_cells = active_cells;

firing_rate_norm = f_normalize(firing_rate, normalize1);

[num_cells, ~] = size(firing_rate_norm);

%% dim reduction with SVD to calulate components number

% [~,S,~] = svd(firing_rate_norm);
% sing_val_sq = diag(S'*S);
% d_explained = sing_val_sq/sum(sing_val_sq)*100;

[d_coeff,~,~,~,d_explained,~] = pca(firing_rate_norm');

dimensionality_total_norm = sum(cumsum(d_explained)<(total_dim_thresh*100));
%figure; plot(d_explained)

%% repeat with not norm
[~,S2,~] = svd(firing_rate);
sing_val_sq2 = diag(S2'*S2);
d_explained2 = sing_val_sq2/sum(sing_val_sq2)*100;
%figure; plot(d_explained)
dimensionality_total = sum(cumsum(d_explained2)<(total_dim_thresh*100));

%% shuff and PCA
max_lamb_shuff = zeros(dim_est_num_reps,1);
dim_total_shuff = zeros(dim_est_num_reps,1);
for n_rep = 1:dim_est_num_reps
    firing_rate_shuff = f_shuffle_data(firing_rate_norm, shuffle_method);
%     [~,s_S,~] = svd(firing_rate_shuff);
%     s_sing_val_sq = diag(s_S'*s_S);
%     s_explained = s_sing_val_sq/sum(s_sing_val_sq)*100;
%     
    [s_coeff,~,~,~,s_explained,~] = pca(firing_rate_shuff');
    
    dim_total_shuff(n_rep) = sum(cumsum(s_explained)<(total_dim_thresh*100));
    max_lamb_shuff(n_rep) = max(s_explained);
end
dimensionality_total_norm_shuff = mean(dim_total_shuff);
% eigenvalues below lower bound plus above upper should
% theoretically equal total number of neurons in all ensembles
%dimensionality_corr = sum(d_explained>prctile(max_lamb_shuff, corr_comp_thresh*100));
%dimensionality_corr = mean(sum(d_explained>max_lamb_shuff'));
%dimensionality_corr = mean(sum(d_explained>max_lamb_shuff'))+std(num_comps);

dimensionality_corr = mean(sum(d_explained>max_lamb_shuff'));
num_comps = ceil(dimensionality_corr);

data_dim_est.dimensionality_total = dimensionality_total;
data_dim_est.dimensionality_first_comp_size = d_explained2(1);
data_dim_est.dimensionality_total_norm = dimensionality_total_norm;
data_dim_est.dimensionality_total_norm_shuff = dimensionality_total_norm_shuff;
data_dim_est.dimensionality_corr = dimensionality_corr;
data_dim_est.num_comps = num_comps;
data_dim_est.d_explained = d_explained(1:num_comps);
%data_dim_est.corr_comp_thresh = corr_comp_thresh;
data_dim_est.num_cells = num_cells;

%%
if plot_stuff
    SI_firing_rate = similarity_index(firing_rate_norm, firing_rate_norm);
    SI_firing_rate_shuff = similarity_index(firing_rate_shuff, firing_rate_shuff);
    
    figure; 
    ax1 = subplot(3,1,1:2);imagesc(firing_rate_norm); axis tight;
    title('Firing rates raster');
    ax2 = subplot(3,1,3);plot(sum(firing_rate_norm)); axis tight;
    linkaxes([ax1,ax2],'x');

    figure; 
    ax1 = subplot(3,1,1:2);imagesc(firing_rate_shuff); axis tight;
    title(['Shuffled rates raster (' shuffle_method ')']);
    ax2 = subplot(3,1,3);plot(sum(firing_rate_shuff)); axis tight;
    linkaxes([ax1,ax2],'x');
    
    
    figure; imagesc(SI_firing_rate);
    title('cell cell similarity');

    figure; imagesc(SI_firing_rate_shuff);
    title('Shuffled cell cell similarity');
    
    figure; hold on;
    ecdf(SI_firing_rate(:));
    ecdf(SI_firing_rate_shuff(:));
    legend('Data', 'Shuffled');
    title('ECDF cell-cell SI');
end


if plot_stuff
    
    figure; 
    ax1 = subplot(2,1,1); hold on; axis tight
    bins = floor(min(d_explained)):0.05:ceil(max(d_explained));
    histogram(d_explained, 'BinEdges', bins);
    histogram(s_explained, 'BinEdges', bins);
    lambda_thresh = mean(max_lamb_shuff);
    line([lambda_thresh lambda_thresh], [ax1.YLim], 'Color', 'r', 'LineStyle', '--')
    title('variance explained dist');
    legend('data', 'shuff', 'thresh');
    ylabel('component number')
    ax2 = subplot(2,1,2); hold on; axis tight
    [f_d, x_d] = ecdf(d_explained);
    [f_s, x_s] = ecdf(s_explained);
    plot(x_d, f_d, 'LineWidth', 2);
    plot(x_s, f_s, 'LineWidth', 2);
    line([lambda_thresh lambda_thresh], [0 1], 'color', 'red', 'LineStyle', '--')
    title('ECDF of variance explained')
    legend('data', 'shuff', 'thresh');
    linkaxes([ax1,ax2],'x');
    ax2.XLim(2) = max(x_d);
    xlabel('Variance per component');
    ylabel('component fraction')
    
    figure; hold on;
    for n_comp = 1:num_comps
        comp_corr = d_coeff(:,n_comp);  
        scatter(n_comp*ones(num_cells,1), comp_corr, '.')
    end
    title('PCA coefficients')
    
    figure; hold on;
    for n_comp = 1:num_comps
        comp_corr = s_coeff(:,n_comp);  
        scatter(n_comp*ones(num_cells,1), comp_corr, '.')
    end
    title('PCA coefficients shuff')
    
    figure;
    ax1 = subplot(3,1,1); hold on;
    plot(d_explained);
    line([0 num_cells], [lambda_thresh, lambda_thresh],'Color','red','LineStyle','--')
    ylabel('Eigenvalues');
    ax2 = subplot(3,1,2:3);
    imagesc(d_coeff);
    ylabel('neuron');
    xlabel('Component')
    linkaxes([ax1,ax2],'x'); axis tight;
    
    for n_pc = 1:3
        figure;
        ax1 = subplot(3,1,1); hold on;
        plot(d_coeff(:,n_pc));
        title(['Principal component ' num2str(n_pc)]);
        ax2 = subplot(3,1,2:3);
        imagesc(d_coeff(:,n_pc)*d_coeff(:,n_pc)');
        title(['PC ' num2str(n_pc) ' outer product']);
        ylabel('neuron'); xlabel('neuron');
        linkaxes([ax1,ax2],'x'); axis tight;
    end
    
end


end