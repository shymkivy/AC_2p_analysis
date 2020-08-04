function data_dim_est = f_ensemble_comp_data_dim2(firing_rate, normalize)
% parameters
shuffle_method = 'scramble'; % 'circ_shift' or 'scramble'
plot_stuff = 0;
var_thresh_prc = .95; % circular shift thresh (95 or 99; from Detecting cell assemblies in large neuronal populations)
%num_comp = 100;
%ensamble_method = 'tca'; % 'PCA', 'AV', 'ICA', 'NMF', 'SPCA', 'tca', 'fa', 'gpfa'
%example_plot = 20;
total_dim_thresh = .7;
%%

if ndims(firing_rate) == 3
    [num_cells, ~, num_trials] = size(firing_rate);
    firing_rate = reshape(firing_rate, num_cells,[]);
else
    [num_cells, num_trials] = size(firing_rate);
end

active_cells = sum(firing_rate,2) > 0;
firing_rate(~active_cells,:) = [];

if normalize
    firing_rate_norm = firing_rate - mean(firing_rate,2);
    firing_rate_norm = firing_rate_norm./std(firing_rate_norm,[],2); 
    %firing_rate_cont(isnan(firing_rate_cont)) = 0;
else
    firing_rate_norm = firing_rate;
end


%% dim reduction with PCA to calulate components number

[~,S,~] = svd(firing_rate_norm);
sing_val_sq = diag(S'*S);
d_explained = sing_val_sq/sum(sing_val_sq)*100;
%figure; plot(d_explained)
dimensionality_total = sum(cumsum(d_explained)<(total_dim_thresh*100));
%[coeff,score,~,~,d_explained,~] = pca(firing_rate_norm');


%% shuff and PCA

num_reps = 20;
data_thresh = zeros(num_reps,1);
dim_total_shuff = zeros(num_reps,1);
for n_rep = 1:num_reps
    firing_rate_shuff = f_shuffle_data(firing_rate_norm, shuffle_method);
    [~,s_S,~] = svd(firing_rate_shuff);
    s_sing_val_sq = diag(s_S'*s_S);
    s_explained = s_sing_val_sq/sum(s_sing_val_sq)*100;
    
    dim_total_shuff(n_rep) = sum(cumsum(s_explained)<(total_dim_thresh*100));
    
    %[~,~,~,~,s_explained,~] = pca(firing_rate_shuff');
    data_thresh(n_rep) = prctile(s_explained, var_thresh_prc*100); % ss_explained or s_explained
end
dimensionality_total_shuff = mean(dim_total_shuff);
% eigenvalues below lower bound plus above upper should
% theoretically equal total number of neurons in all ensembles
dimensionality_corr = sum(sum(d_explained>data_thresh'))/num_reps;

num_comps = max(ceil(dimensionality_total),1);
data_dim_est.dimensionality_total = dimensionality_total;
data_dim_est.dimensionality_total_shuff = dimensionality_total_shuff;
data_dim_est.dimensionality_corr = dimensionality_corr;
data_dim_est.num_comps = num_comps;
data_dim_est.d_explained = d_explained(1:num_comps);
data_dim_est.var_thresh_prc = var_thresh_prc;
data_dim_est.num_cells = num_cells;
data_dim_est.num_trials = num_trials;

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
    bins = floor(min(d_explained)):0.05:ceil(max(d_explained));
    figure; hold on;
    histogram(d_explained, 'BinEdges', bins);
    histogram(s_explained, 'BinEdges', bins);
    title('variance explained dist');
    legend('data', 'shuff');

    figure; hold on;
    ecdf(d_explained);
    ecdf(s_explained);
    line([0 4], [var_thresh_prc var_thresh_prc], 'color', 'red')
    title('ECDF of variance explained')
    legend('data', 'shuff', [num2str(var_thresh_prc*100) '% thresh']);
    
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
    line([0 num_cells], [data_thresh, data_thresh],'Color','red','LineStyle','--')
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