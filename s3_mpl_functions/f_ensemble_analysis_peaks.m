function data_dim_est = f_ensemble_analysis_peaks(trial_data,trial_types, params, ops)
% parameters
cond_name = params.cond_name;
n_dset = params.n_dset;

normalize_firing_rates = false;
shuffle_method = 'scramble'; % 'circ_shift' or 'scramble'
plot_stuff = 0;
plot_stuff_extra = 0;
plot_ens_details = 0;
var_thresh_prc = ops.ensemb.pca_var_thresh; % circular shift thresh (95 or 99; from Detecting cell assemblies in large neuronal populations)
ensamble_method = ops.ensemb.method; % 'PCA', 'AV', 'ICA', 'NMF', 'SPCA', 'tca', 'fa', 'gpfa'
% so far NMF and ICA seems to work

%%

num_trials = numel(trial_types);
if ndims(trial_data) == 3
    [num_cells, ~, num_trials] = size(trial_data);
    trial_data = reshape(trial_data, num_cells,[]);
end  
num_bins = size(trial_data,2)/num_trials;

active_cells = sum(trial_data,2) > 0;
trial_data(~active_cells,:) = [];

if normalize_firing_rates
    firing_rate_norm = trial_data - mean(trial_data,2);
    firing_rate_norm = firing_rate_norm./std(firing_rate_norm,[],2); 
    %firing_rate_cont(isnan(firing_rate_cont)) = 0;
else
    firing_rate_norm = trial_data;
end
[num_cells, ~] = size(firing_rate_norm);
%firing_rate_norm_3d = reshape(firing_rate_norm, num_cells, num_bins, num_trials);

%% dim reduction with PCA to calulate components number

[d_coeff,d_score,~,~,d_explained,d_mu] = pca(firing_rate_norm');


%% shuff and PCA
num_reps = 100;
data_thresh = zeros(num_reps,1);
for n_rep = 1:num_reps
    firing_rate_shuff = f_shuffle_data(firing_rate_norm, shuffle_method);
    [~,~,~,~,s_explained,~] = pca(firing_rate_shuff');
    data_thresh(n_rep) = prctile(s_explained, var_thresh_prc*100); % ss_explained or s_explained
end

% eigenvalues below lower bound plus above upper should
% theoretically equal total number of neurons in all ensembles

dimensionality = sum(sum(d_explained>data_thresh'))/num_reps;
num_comps = max(ceil(dimensionality),1);

data_dim_est.dimensionality = dimensionality;
data_dim_est.num_comps = num_comps;
data_dim_est.d_explained = d_explained(num_comps);
data_dim_est.var_thresh_prc = var_thresh_prc;
data_dim_est.num_cells = num_cells;
data_dim_est.num_trials = num_trials;

%%
SI_firing_rate = similarity_index(firing_rate_norm, firing_rate_norm);
SI_firing_rate_shuff = similarity_index(firing_rate_shuff, firing_rate_shuff);

firing_rate_LR = d_coeff(:,1:num_comps)*d_score(:,1:num_comps)' + d_mu';
SI_firing_rate_LR = similarity_index(firing_rate_LR, firing_rate_LR);

%%
%use_projection_mat = 1;
% if use_projection_mat
%     firing_rate_norm_dred = 
% end

%% real data 

ensamble_method = 'svd';

num_ens_comps = num_comps;
if strcmpi(ensamble_method, 'nmf')
    num_ens_comps = round(num_comps*1.5);
end

[dred_factors1, ~] = f_dred_train2(firing_rate_LR, num_ens_comps, ensamble_method, 0);
[coeffs, scores] = f_dred_get_coeffs(dred_factors1);

%% Visualize traces
tn_to_dred = params.tn_to_dred;
tt_to_dred = params.tt_to_dred;

X = scores';

num_plots = ceil(num_comps/3);
for n_plt = 1:num_plots
    figure; hold on;
    lg1 = cell(numel(tt_to_dred),1);
    dims1 = rem([1 2 3]+(n_plt-1)*3-1, num_comps)+1;
    for n_tt = 1:numel(tt_to_dred)
        plot3(X(trial_types == tt_to_dred(n_tt),dims1(1)), X(trial_types == tt_to_dred(n_tt),dims1(2)),X(trial_types == tt_to_dred(n_tt),dims1(3)), 'Color', ops.colors_list{n_tt},'Marker', 'o', 'LineStyle', 'none'); 
        axis tight;
        lg1{n_tt} = ops.context_types_labels{ops.context_types_all==tt_to_dred(n_tt)};
    end
    xlabel(sprintf('pc %d', dims1(1)));
    ylabel(sprintf('pc %d', dims1(2)));
    zlabel(sprintf('pc %d', dims1(3)));
    title(sprintf('%s, dset%d, trial types pcs scores(trials)',params.cond_name, params.n_dset));
    legend(lg1)
    % %quiver3 for gradient plotting
end


params2.method = 'ward'; % cosine, ward
params2.metric = 'euclidean'; % cosine squaredeuclidean
params2.plot_sm = 0;
hclust_out = f_hcluster_trial3(X, params2);
figure;
gscatter(X(:,1),X(:,2),hclust_out.clust_ident);


params3.metric = 'sqEuclidean'; % cosine squaredeuclidean
params3.RegularizationValue = 0.01;
gmmclust_out = f_gmmcluster_trial(X, params3);
figure;
gscatter(X(:,1),X(:,2),gmmclust_out.clust_ident);


num_plots = ceil(num_comps/3);
for n_plt = 1:num_plots
    figure; hold on;
    dims1 = rem([1 2 3]+(n_plt-1)*3-1, num_comps)+1;
    for n_tt = 1:numel(tt_to_dred)
        plot3(coeffs(:,dims1(1)), coeffs(:,dims1(2)),coeffs(:,dims1(3)), 'Color', ops.colors_list{n_tt},'Marker', 'o', 'LineStyle', 'none'); 
        axis tight;
    end
    xlabel(sprintf('pc %d', dims1(1)));
    ylabel(sprintf('pc %d', dims1(2)));
    zlabel(sprintf('pc %d', dims1(3)));
    title(sprintf('%s, dset%d, cell types pcs coeffs(cells)',params.cond_name, params.n_dset));
    % %quiver3 for gradient plotting
end



% figure; imagesc(firing_rate_norm)
% for n_comp = 1:num_comps
%     figure;
%     plot(coeffs(:,n_comp));
%     xlabel('cells')
%     title(sprintf('Coeff %d', n_comp));
%     figure;
%     plot(scores(n_comp,:));
%     title(sprintf('Score %d', n_comp));
%     xlabel('trials')
% end


num_rep = 10000;
figure; hold on;
for n_tt = 1:numel(tt_to_dred)
    tt_list = find(trial_types == tt_to_dred(n_tt));
    dist1 = zeros((numel(tt_list)-1),1);
    for n_rep = 1:(numel(tt_list)-1)
        dist1(n_rep) = (scores(:,tt_list(n_rep))/norm(scores(:,tt_list(n_rep))))'*(scores(:,tt_list(n_rep+1))/norm(scores(:,tt_list(n_rep+1))));
    end
    [f, x] = ecdf(dist1);
    plot(x, f, 'color', ops.colors_list{n_tt}, 'lineWidth', 2);
end
for n_tt = 1:numel(tt_to_dred)
    dist1 = zeros(num_rep,3);
    tt_list = find(trial_types == tt_to_dred(n_tt));
    for n_rep = 1:num_rep
        samp_tr = randsample(tt_list,2);
        dist1(n_rep,n_tt) = (scores(:,samp_tr(1))/norm(scores(:,samp_tr(1))))'*(scores(:,samp_tr(2))/norm(scores(:,samp_tr(2))));
    end
    [f, x] = ecdf(dist1(:,n_tt));
    plot(x, f, 'color', ops.colors_list{n_tt}, 'LineStyle', '--', 'lineWidth', 2);
end
title(sprintf('%s, dset%d, Onset population trajectory structure ECDF',params.cond_name, params.n_dset));
legend('cont', 'red', 'dev', 'rand', 'Location', 'northwest');
xlabel('cosine similarity betweeen trials');
ylabel('Fraction');


if strcmpi(ensamble_method, 'tca')
    tn_seq = sum((trial_types == tt_to_dred').*tn_to_dred,2);
    [~, tt_idx] = sort(tn_seq, 'ascend');
    trial_types_dred_sort = trial_types(tt_idx);

    %f_viz_ktensor_trials(dred_factors1.dred_factors.t_factors, trial_types_dred_sort)
    %f_dred_plot_factors(dred_data_list3,trial_types_dred, test_data_ind);

    t_factors_sort = dred_factors1.dred_factors.t_factors;%.U{3}
    t_factors_sort.U{3} = t_factors_sort.U{3}(tt_idx,:);

    %if ops.dred_params.dred_mmn
    c_map = ops.plot_params.c_map_mmn;
    %c_map = [.6 .6 .6; .8 .8 .8; .5 .5 1; 1 .5 .5; .7 .7 1; 1 .7 .7];
    %else
        %c_map = jet(numel(tt_to_dred));
    %end
    trial_types_dred_sort_colors = zeros(numel(trial_types_dred_sort),3);
    for n_tt = 1:numel(tt_to_dred)
        trial_types_dred_sort_colors(trial_types_dred_sort == tt_to_dred(n_tt),:) = trial_types_dred_sort_colors(trial_types_dred_sort == tt_to_dred(n_tt),:)+c_map(n_tt,:);
    end
    trial_image_sort = reshape(trial_types_dred_sort_colors,1,numel(trial_types_dred_sort),3);
    %figure; imagesc(trial_image_sort);

    trial_types_dred_reg_colors = zeros(numel(trial_types),3);
    for n_tt = 1:numel(tt_to_dred)
        trial_types_dred_reg_colors(trial_types == tt_to_dred(n_tt),:) = trial_types_dred_reg_colors(trial_types == tt_to_dred(n_tt),:)+c_map(n_tt,:);
    end
    trial_image_reg = reshape(trial_types_dred_reg_colors,1,numel(trial_types),3);
    %figure; imagesc(trial_image_reg);

    %
    %viz_ktensor_mmn(dred_factors1.dred_factors.t_factors,trial_image_reg)%, 'Sameylims', false(num_ens_comps,1)
    %suptitle([cond_name ' dset' num2str(n_dset) ' unsorted tca']);
    viz_ktensor_mmn(t_factors_sort,trial_image_sort) %, 'Sameylims', false(num_ens_comps,1)
    suptitle([cond_name ' dset' num2str(n_dset) ' sorted tca']);

    %figure; plot(dred_factors1.dred_factors.t_factors.lambda)
    %title([cond_name ' dset' num2str(n_dset) ' TCA lambda']);
end


%% use shuffle data to get thresholds for ensemble acception


shuff_rep = 10;
shuff_data = cell(shuff_rep,1);
for n_rep = 1:shuff_rep
    firing_rate_LR_shuff = f_shuffle_data(firing_rate_LR, 'circ_shift');
    train_done = 0;
    while ~train_done
        try
            [dred_factors_shuff, ~] = f_dred_train2(firing_rate_LR_shuff, num_ens_comps, ensamble_method, 0);
            train_done = 1;
        catch
            disp('Error train, will repeat');
        end
    end
    coeffs_shuff = f_dred_get_coeffs(dred_factors_shuff);
    shuff_data{n_rep} = coeffs_shuff;
end
shuff_data1 = cat(1,shuff_data{:});
if ~sum(shuff_data1(:)<0)
    ens_thresh = prctile(shuff_data1(:), 95);
else
    ens_thresh = prctile(shuff_data1(:), [2.5 97.5]);
end


%% extract ensambles
ens = f_extract_ensamble_from_scores(dred_factors1, ens_thresh);
%suptitle([cond_name ' dset' num2str(n_dset) ' ' ensamble_method]);
%ens = f_extract_ensamble_from_scores(dred_factors_shuff);
data_dim_est.ensambles = ens;

%[~, e_scores] = f_ens_get_coeffs(dred_factors1);


%% plot 

% thresholds
if plot_ens_details
    figure; hold on;
    for n_ens = 1:num_ens_comps
        subplot(5, ceil(num_ens_comps/5),n_ens); hold on;
        for n_comp = 1:num_comps
            [f, x] = ecdf(shuff_data1(:,n_comp));%ksdensity(shuff_data{1}(:,n_comp), 'Bandwidth', 0.1);
            plot(x, f, 'Color', [.6 .6 .6], 'LineWidth', 0.1)
        end
        [f, x] = ecdf(shuff_data1(:));
        ax1 = plot(x, f, 'r', 'LineWidth', 2);
        [f, x] = ecdf(coeffs(:,n_ens));
        ax2 = plot(x, f, 'k', 'LineWidth', 2);
        for n_thr = 1:numel(ens_thresh)
            ax3 = line([ens_thresh(n_thr) ens_thresh(n_thr)], [0 1], 'Color', 'r','LineStyle','--');
        end
        axis tight;
        if n_ens == 1
            legend([ax1 ax2 ax3], {'Shuffle', 'Ensamble', 'ens thresh'});
        end
        title(num2str(n_ens))
    end
end

% ensambles
if plot_ens_details
    % sort according to context and make color map for trial labels
    if numel(intersect(params.ctx_mmn,trial_types))>= 6
        c_map_dd = [.6 .6 .6; .6 .6 1; 1 .6 .6; .8 .8 .8; .6 .6 1; 1 .6 .6];
        % extract ctx trials
        trial_types_ctx_index = logical(sum(trial_types == params.ctx_mmn',2));
        trial_types_ctx = trial_types(trial_types_ctx_index);
        firing_rate_norm_3d_ctx = firing_rate_norm_3d(:,:,trial_types_ctx_index);
        
        [~, trial_types_idx_sort] = sort(sum((trial_types_ctx == params.ctx_mmn').*(1:6),2), 'ascend');
        trial_types_ctx_sort = trial_types_ctx(trial_types_idx_sort);
        firing_rate_norm_3d_ctx_sort = reshape(firing_rate_norm_3d_ctx(:,:,trial_types_idx_sort),num_cells, []);

        trace_map_ctx = if_generate_trial_trace_map(trial_types_ctx, params.ctx_mmn', num_bins, c_map_dd);
        trace_map_ctx_sort = if_generate_trial_trace_map(trial_types_ctx_sort, params.ctx_mmn', num_bins, c_map_dd);
    end
    
    
    cont_trials = intersect(1:10,trial_types);
    if numel(cont_trials)
        c_map_cont = jet(numel(cont_trials));
        trial_types_cont_index = logical(sum(trial_types == 1:10,2));
        trial_types_cont = trial_types(trial_types_cont_index);
        firing_rate_norm_3d_cont = firing_rate_norm_3d(:,:,trial_types_cont_index);
        
        [~, trial_types_idx_sort] = sort(sum((trial_types_cont == 1:10).*(1:10),2), 'ascend');
        trial_types_sort_cont = trial_types_cont(trial_types_idx_sort);
        firing_rate_norm_3d_cont_sort = reshape(firing_rate_norm_3d(:,:,trial_types_idx_sort),num_cells, []);

        trace_map_cont = if_generate_trial_trace_map(trial_types_cont, cont_trials', num_bins, c_map_cont);
        trace_map_cont_sort = if_generate_trial_trace_map(trial_types_sort_cont, cont_trials', num_bins, c_map_cont);
    end
    
    for n_ens = 1:numel(ens)
        coeffs_ens = zeros(num_cells,1);
        coeffs_ens(ens(n_ens).ens_cell_num) = ens(n_ens).cell_comp_mag;
        ens_traces = coeffs_ens*ens(n_ens).ens_score;
        if plot_stuff_extra
            figure;
            ax1 = subplot(4,1,1);
            imagesc(firing_rate_norm(ens(n_ens).ens_cell_num,:));
            title(['Full rank data; ensamble ' num2str(n_ens)]);
            ax2 = subplot(5,1,2);
            plot(mean(firing_rate_norm(ens(n_ens).ens_cell_num,:)));
            title('mean full rank')
            ax3 = subplot(5,1,3);
            imagesc(firing_rate_LR(ens(n_ens).ens_cell_num,:));
            title('Low rank data');
            ax4 = subplot(5,1,4);
            imagesc(ens_traces(ens(n_ens).ens_cell_num,:));
            title('Reconstructed from sig');
            ax5 = subplot(5,1,5);
            plot(ens(n_ens).ens_score); 
            title([ensamble_method ' ensamble']);
            linkaxes([ax1,ax2, ax3, ax4, ax5],'x');
            axis tight;
        end
        
        % plot control data
        if numel(cont_trials)>2
            figure;
            [m,n] = f_give_subplotdims(numel(cont_trials));
            sp = cell(numel(cont_trials),1);
            for n_tr = 1:numel(cont_trials)
                sp{n_tr} = subplot(m+3,n,[n_tr]);
            end
            ens_data = reshape(trial_data(ens(n_ens).ens_cell_num,:),numel(ens(n_ens).ens_cell_num), [], num_trials);
            trial_ave = f_mpl_trial_average(ens_data,trial_types, cont_trials, 'none');
            f_plot_trials(trial_ave, cont_trials, params.trial_t, sp)
            sp1 = subplot(2*(m+3),n,2*m*n+(1:2*n));
            h = imagesc(firing_rate_norm_3d_cont(ens(n_ens).ens_cell_num,:));
            if_add_trials_to_plot(sp1, h, trace_map_cont)
            title(['cont trials ensamble ' num2str(n_ens)]);
            sp2 = subplot(2*(m+3),n,2*m*n+2*n+[1 n]);
            plot(mean(firing_rate_norm_3d_cont(ens(n_ens).ens_cell_num,:)));
            linkaxes([sp1,sp2],'x'); axis tight;
            sp3 = subplot(2*(m+3),n,2*m*n+3*n+(1:2*n));
            h = imagesc(firing_rate_norm_3d_cont_sort(ens(n_ens).ens_cell_num,:));
            if_add_trials_to_plot(sp3, h, trace_map_cont_sort)
            title('Sorted ctx trials')
            sp4 = subplot(2*(m+3),n,2*m*n+5*n+[1 n]);
            plot(mean(firing_rate_norm_3d_cont_sort(ens(n_ens).ens_cell_num,:)));
            linkaxes([sp3, sp4],'x'); axis tight;
            
        end

        % plot mmn data
        if numel(intersect(params.ctx_mmn,trial_types))>= 6
            figure;
            sp = cell(2,1);
            sp{1} = subplot(8,2,[1 3]);
            sp{2} = subplot(8,2,[2 4]);
            ens_data = reshape(trial_data(ens(n_ens).ens_cell_num,:),numel(ens(n_ens).ens_cell_num), [], num_trials);
            trial_ave = f_mpl_trial_average(ens_data,trial_types, params.ctx_mmn, 'none');
            f_mpl_plot_dd(trial_ave, ones(size(ens_data,1),numel(params.ctx_mmn)), params.trial_t, ops,sp)
            sp1 = subplot(8,2,5:8);
            h = imagesc(firing_rate_norm_3d_ctx(ens(n_ens).ens_cell_num,:));
            if_add_trials_to_plot(sp1, h, trace_map_ctx)
            title(['ctx trials ensamble ' num2str(n_ens)])
            sp2 = subplot(8,2,[9 10]);
            plot(mean(firing_rate_norm_3d_ctx(ens(n_ens).ens_cell_num,:)));
            linkaxes([sp1,sp2],'x'); axis tight;
            sp3 = subplot(8,2,11:14);
            h = imagesc(firing_rate_norm_3d_ctx_sort(ens(n_ens).ens_cell_num,:));
            if_add_trials_to_plot(sp3, h, trace_map_ctx_sort)
            title('Sorted ctx trials')
            sp4 = subplot(8,2,[15 16]);
            plot(mean(firing_rate_norm_3d_ctx_sort(ens(n_ens).ens_cell_num,:)));
            linkaxes([sp3,sp4],'x'); axis tight;
        end
    end
end



%% need to edit from here
% 
% for num_pc = 1:example_plot
%     n_pc = plot_cells(num_pc);
%     for ii = 1:2
%         if ~isempty(ens_cells{n_pc}{ii})
%             figure;
%             subplot(4,1,1:2);
%             imagesc(firing_rate_norm(ens_cells{n_pc}{ii},:));
%             title(['PC ' num2str(n_pc)])
%             subplot(4,1,3);
%             plot(sum(firing_rate_norm(ens_cells{n_pc}{ii},:)));
%             subplot(4,1,4);
%             plot(ens_act(n_pc,:))
%         end
%     end
% end
% 
% 
% for n_pc = 1:10
%     figure; 
%     ax1 = subplot(4,5,1:4);
%     plot(proj_dist);
%     title(['PC ' num2str(n_pc) ' trace']);
%     ax2 = subplot(4,5,[6:9,11:14]);
%     imagesc(firing_rate_norm);
%     title('Firing raster');
%     ylabel('neurons');
%     ax3 = subplot(4,5,[10,15]);
%     plot(w,1:numel(w));
%     ax3.YDir = 'reverse'; axis tight;
%     title(['PC ' num2str(n_pc) ' neurons']);
%     ax4 = subplot(4,5,16:19);
%     plot(sum(firing_rate_norm));
%     title('Population firing rate');
%     xlabel('time');
%     linkaxes([ax1,ax2, ax4],'x'); axis tight;
%     linkaxes([ax2,ax3],'y'); 
% end
% 
% 
% for n_pc = 1:50
%     if sum(A(:,n_pc)>median(A(:,n_pc))+3*std(A(:,n_pc)))
%         disp(['PC ' num2str(n_pc) ', ' num2str(sum(A(:,n_pc)>median(A(:,n_pc))+2*std(A(:,n_pc)))) ' positive components'])
%     end
%     if sum(A(:,n_pc)<median(A(:,n_pc))-3*std(A(:,n_pc)))
%         disp(['PC ' num2str(n_pc) ', ' num2str(sum(A(:,n_pc)<median(A(:,n_pc))-2*std(A(:,n_pc)))) ' negative components'])
%     end
% end
% 
% 
% n_pc = 6;
% figure; 
% subplot(3,1,1:2);
% imagesc(firing_rate_norm(A(:,n_pc)<median(A(:,n_pc))-2*std(A(:,n_pc)),:))
% subplot(3,1,3);
% plot(sum(firing_rate_norm(A(:,n_pc)<median(A(:,n_pc))-2*std(A(:,n_pc)),:)))
% 
% proj2 = A*icasig;
% 
% figure; hold on;
% plot(proj(50,:))
% plot(proj2(50,:))
% 
% 
% figure;
% plot(Out1(20,:)');
% title('FastICA components')
% 
% figure; plot(Out2(:,15))
% 
% figure; plot(Out3(1,:))
% 
% 
% figure; plot(M(:,3))
% 
% figure; histogram(M(:,1))
% figure; histogram(M(:,12))
% figure; plot(M(:,12))
% 
% figure;
% scatter3(M(:,1),M(:,2),M(:,3), 'o')
% 
% figure; imagesc(corr_mat); title('corr mat');
% figure; imagesc(N); title('Corr projected on PC subspace');
% figure; imagesc(M); title('interaction matrix');
% 
% disp([num2str(num_comps) ' PCA components above ' num2str(var_thresh_prc*100) '% thresh']);
% disp([num2str(sum(isnan(SI_firing_rate(:)))) ' SI values are NaN']);
% %SI_firing_rate(isnan(SI_firing_rate)) = 0;
% 
% [d_U,d_S,d_V] = svd(SI_firing_rate);
% [s_U,s_S,s_V] = svd(SI_firing_rate_shuff);
% 
% 
% figure; hold on;
% histogram(d_U(:,1))
% histogram(d_U(:,2))
% histogram(d_U(:,3))
% histogram(d_U(:,4))
% histogram(d_U(:,5))
% histogram(d_U(:,6))
% histogram(d_U(:,7))
% histogram(d_U(:,8))
% histogram(d_U(:,9))
% title('SVD coefficients distributions')
% 
% figure; hold on;
% plot(std(s_V))
% plot(std(s_V'))
% plot(std(d_V))
% plot(std(d_V'))
% title('STD of population eigenvectors');
% legend('shuffle V', 'shuffle VT', 'V', 'VT');
% 
% 
% 
% % error needs to be positive definite
% %[d_lambda,d_psi,d_T,d_stats] = factoran(firing_rate_cont,10);
% %[s_lambda,s_psi,s_T,s_stats] = factoran(firing_rate_cont_shuffcirc,10);
% 
% 
% 
% Mdl = rica(firing_rate_norm,num_comp);
% 
% figure; hold on;
% plot(Mdl.TransformWeights(:,plot_comp_num))
% title('RICA components')
% 
% figure; hold on;
% for n_comp = 1:num_comp
%     comp_corr = firing_rate_norm*Mdl.TransformWeights(:,n_comp);  
%     scatter(n_comp*ones(num_cells,1), comp_corr, '.')
% end
% title('RICA coefficients')
% 
% 
% [Out1, Out2, Out3] = fastica(firing_rate_norm,'lastEig', num_comp ); % 'lastEig' 
% figure;
% plot(Out1(plot_comp_num,:)');
% title('FastICA components')
% 
% 
% 
% 
% 
% figure; hold on;
% plot(mean(firing_rate_norm(comp_corr<3*std(comp_corr),:)))
% plot(Mdl.TransformWeights(:,comp))
% 
% 
% figure; histogram(x)
% 
% 
% x = d_U(:,1:10)'*SI_firing_rate;
% figure; plot(x(1,:))
% figure; plot(x(2,:))
% 
% 
% figure; histogram(diag(s_S))
% 
% cov_mat1 = firing_rate_norm*firing_rate_norm';
% 
% figure; imagesc(cov_mat1)
% figure; imagesc(SI_firing_rate)
% 
% x = c_U*c_U
% 
% [c_U,c_S,c_V] = svd(cov_mat1);
% 
% figure; plot(c_U(:,1)'*c_U(:,:))
% 
%% plots 
if plot_stuff
    figure; 
    ax1 = subplot(3,1,1:2);imagesc(firing_rate_norm);
    caxis([0 20]);
    title(['Firing rates raster ' cond_name '; dset ' num2str(n_dset)]);
    ax2 = subplot(3,1,3);plot(sum(firing_rate_norm));
    linkaxes([ax1,ax2],'x');
    axis tight;
    
    
    figure; 
    ax3 = subplot(3,1,1:2);imagesc(firing_rate_shuff);
    title(['Shuffled rates raster (' shuffle_method ') ' cond_name '; dset ' num2str(n_dset)]);
    ax4 = subplot(3,1,3);plot(sum(firing_rate_shuff));
    linkaxes([ax3,ax4],'x');
    axis tight;
    y_min = min([ax2.YLim ax4.YLim]);
    y_max = max([ax2.YLim ax4.YLim]);
    ax2.YLim = [y_min y_max];
    ax4.YLim = [y_min y_max];
    
    
    figure; 
    imagesc(SI_firing_rate);  %  - eye(size(SI_firing_rate,1))
    title(['cell cell similarity ' cond_name '; dset ' num2str(n_dset)]);
    caxis([-.2 .8]);
    axis equal tight;
    xlabel('Cells');
    ylabel('Cells');
    colorbar;
    
    %Z = pdist(firing_rate_norm, 'cosine'); %linkage(Z,'ward')
    figure;
    [~, ~, dend_order] = dendrogram(linkage(firing_rate_norm,'ward'), 1000,'ColorThreshold',0.65);
    Z = pdist(firing_rate_norm(dend_order,:), 'cosine');
    imagesc(1-squareform(Z));
    axis image;
    title(['Cell to cell similarity sorted ward ' cond_name '; dset ' num2str(n_dset)]);
    caxis([-.2 .8]);
    axis equal tight;
    xlabel('Cells');
    ylabel('Cells');
    colorbar

    figure;
    imagesc(SI_firing_rate_shuff);
    title(['Shuffled cell cell similarity ' cond_name '; dset ' num2str(n_dset)]);
    caxis([-.2 .8]);
    axis equal tight;
    xlabel('Cells');
    ylabel('Cells');
    colorbar
    
    figure; hold on;
    ecdf(SI_firing_rate(:));
    ecdf(SI_firing_rate_shuff(:));
    legend('Data', 'Shuffled');
    title(['ECDF cell-cell SI ' cond_name '; dset ' num2str(n_dset)]);
end

if plot_stuff
    bins = floor(min(d_explained)):0.1:ceil(max(d_explained));
    figure; hold on;
    histogram(d_explained, 'BinEdges', bins);
    histogram(s_explained, 'BinEdges', bins);
    title(['variance explained dist ' cond_name '; dset ' num2str(n_dset)]);
    legend('data', 'shuff');
    xlabel('eigenvalue magnitude')
    ylabel('count')
    
    figure; hold on;
    [f, x] = ecdf(d_explained);
    plot(x, f, 'LineWidth', 2)
    [f, x] = ecdf(s_explained);
    [~, thr_idx] = min(abs(f - var_thresh_prc));
    x_thresh = x(thr_idx);
    plot(x, f, 'LineWidth', 2)
    line([x_thresh x_thresh], [0 1], 'color', 'red', 'LineStyle', '--')
    title(['ECDF of variance explained ' cond_name '; dset ' num2str(n_dset)])
    legend('data', 'shuff', [num2str(var_thresh_prc*100) '% thresh']);
    xlabel('variance magnitude');
    ylabel('percentile');
    
    
    figure; hold on;
    [f, x] = ksdensity(d_explained);
    plot(x, f, 'LineWidth', 2);
    [f, x] = ksdensity(s_explained);
    plot(x, f, 'LineWidth', 2);
    [~, thr_idx] = min(abs(f - var_thresh_prc));
    line([x_thresh x_thresh], [0 1], 'color', 'red', 'LineStyle', '--')
    title(['Eigenvalue distribution ' cond_name '; dset ' num2str(n_dset)]);
    xlabel('eigenval magnitude');
    ylabel('count');
    legend('data', 'shuff', [num2str(var_thresh_prc*100) '% thresh']);
    
    
    % plot low rank data
    
    figure; 
    ax1 = subplot(3,1,1:2);imagesc(firing_rate_LR);
    caxis([0 20]);
    title(['Firing rates Low Rank ' cond_name '; dset ' num2str(n_dset) '; ' num2str(num_comps) ' comp; ' num2str(round(sum(d_explained(1:num_comps)))) '% var explained']);
    ax2 = subplot(3,1,3);plot(sum(firing_rate_LR));
    linkaxes([ax1,ax2],'x');
    axis tight;

    figure; 
    imagesc(SI_firing_rate_LR);  %  - eye(size(SI_firing_rate,1))
    title(['cell cell similarity Low Rank' cond_name '; dset ' num2str(n_dset)]);    
    axis equal tight;
    caxis([-.2 .8]);
    xlabel('Cells');
    ylabel('Cells');
    
    figure;
    [~, ~, dend_order] = dendrogram(linkage(SI_firing_rate_LR,'ward'), 1000,'ColorThreshold',0.65);
    Z = pdist(firing_rate_norm(dend_order,:), 'cosine');
    imagesc(1-squareform(Z));
    axis image;
    title(['Cell to cell similarity Low Rank sorted ward ' cond_name '; dset ' num2str(n_dset)]);
    caxis([-.2 .8]);
    axis equal tight;
    xlabel('Cells');
    ylabel('Cells');
    colorbar;
end

if plot_stuff_extra
    
    figure; hold on;
    for n_comp = 1:num_comps
        comp_corr = d_coeff(:,n_comp);  
        scatter(n_comp*ones(num_cells,1), comp_corr, '.')
    end
    title(['PCA coefficients ' cond_name '; dset ' num2str(n_dset)])
    
    figure; hold on;
    for n_comp = 1:num_cells
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

%% some clustering
% Z = pdist(firing_rate_norm, 'cosine');
% figure;
% [~, ~, dend_order] = dendrogram(linkage(Z,'average'), 1000,'ColorThreshold',0.65);
% Z = pdist(firing_rate_norm(dend_order,:), 'cosine');
% imagesc(1-squareform(Z));
% axis image;
% title('Cell to cell similarity sorted average');



% figure;
% imagesc(1-squareform(Z));
% title('Cell to cell similarity all ctx');
% xlabel('Cells');
% ylabel('Cells');
% figure;
% dendrogram(linkage(Z,'average'), 1000,'ColorThreshold',0.65);
% title('Cell to cell similarity all ctx');


end

function if_add_trials_to_plot(axis_h, plot_h, trials_trace)

subplot(axis_h)
cdata1 = plot_h.CData;
cmap1 = colormap;
cmin = min(cdata1(:));
cmax = max(cdata1(:));
m = length(cmap1);
index = fix((cdata1-cmin)/(cmax-cmin)*m)+1; 
RGB = ind2rgb(index,cmap1);
imagesc([RGB; trials_trace]);

end

function trace_map = if_generate_trial_trace_map(trial_types_seq, trials_type, num_bins, c_map_dd)

trace_map_temp = sum(repmat(trial_types_seq == trials_type,1,1,3).*reshape(c_map_dd, 1, numel(trials_type), 3),2);
trace_map = reshape(repmat(reshape(trace_map_temp,1,[],3),num_bins,1,1),1,[],3);

end
