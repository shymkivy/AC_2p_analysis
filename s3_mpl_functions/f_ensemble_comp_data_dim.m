function data_dim_est = f_ensemble_comp_data_dim(firing_rate)
% parameters
normalize_firing_rates = true;
shuffle_method = 'circ_shift'; % 'circ_shift' or 'scramble'
plot_stuff = 0;
var_thresh_prc = .95; % circular shift thresh (95 or 99; from Detecting cell assemblies in large neuronal populations)
num_comp = 100;
ensamble_method = 'tca'; % 'PCA', 'AV', 'ICA', 'NMF', 'SPCA', 'tca', 'fa', 'gpfa'
example_plot = 20;

%%

if ndims(firing_rate) == 3
    [num_cells, ~, num_trials] = size(firing_rate);
    firing_rate = reshape(firing_rate, num_cells,[]);
end

active_cells = sum(firing_rate,2) > 0;
firing_rate(~active_cells,:) = [];

if normalize_firing_rates
    firing_rate_norm = firing_rate - mean(firing_rate,2);
    firing_rate_norm = firing_rate_norm./std(firing_rate_norm,[],2); 
    %firing_rate_cont(isnan(firing_rate_cont)) = 0;
else
    firing_rate_norm = firing_rate;
end
[num_cells, num_frames] = size(firing_rate_norm);


%% random shuffle in time
firing_rate_shuff = zeros(num_cells, num_frames);
if strcmp(shuffle_method, 'circ_shift')
    for n_cell = 1:num_cells
        firing_rate_shuff(n_cell,:) = circshift(firing_rate_norm(n_cell,:),ceil(rand(1)*num_frames));
    end
elseif strcmp(shuffle_method, 'scramble')
    % or do also scramble shuffle
    for n_cell = 1:num_cells
        firing_rate_shuff(n_cell,:) = firing_rate_norm(n_cell, randperm(num_frames,num_frames));
    end
end


SI_firing_rate = similarity_index(firing_rate_norm, firing_rate_norm);
SI_firing_rate_shuff = similarity_index(firing_rate_shuff, firing_rate_shuff);


if plot_stuff
    figure; 
    ax1 = subplot(3,1,1:2);imagesc(firing_rate_norm);
    title('Firing rates raster');
    ax2 = subplot(3,1,3);plot(sum(firing_rate_norm));
    linkaxes([ax1,ax2],'x');

    figure; 
    ax1 = subplot(3,1,1:2);imagesc(firing_rate_shuff);
    title(['Shuffled rates raster (' shuffle_method ')']);
    ax2 = subplot(3,1,3);plot(sum(firing_rate_shuff));
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


%% dim reduction with PCA to calulate components number

[d_coeff,d_score,~,~,d_explained,~] = pca(firing_rate_norm');
[s_coeff,s_score,~,~,s_explained,~] = pca(firing_rate_shuff');


% eigenvalues below lower bound plus above upper should
% theoretically equal total number of neurons in all ensembles
data_thresh = prctile(s_explained, var_thresh_prc*100); % ss_explained or s_explained
num_comps = sum(d_explained>data_thresh);

data_dim_est.num_comps = num_comps;
data_dim_est.d_explained = d_explained(logical(1:num_comps));
data_dim_est.var_thresh_prc = var_thresh_prc;
data_dim_est.num_cells = num_cells;
data_dim_est.num_trials = num_trials;


%rec_corr = d_coeff*diag(d_explained)*d_coeff';
%figure; imagesc(rec_corr)

%d_corr = firing_rate_cont*firing_rate_cont';
%figure; imagesc(d_corr)
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
    for n_comp = 1:num_comp
        comp_corr = d_coeff(:,n_comp);  
        scatter(n_comp*ones(num_cells,1), comp_corr, '.')
    end
    title('PCA coefficients')
    
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


%%
% use_projection_mat = 1;
% if use_projection_mat
%     firing_rate_norm_dred = 
% end
%ensamble_method = 'spca';

%[dred_factors1, ~] = f_dred_train2(firing_rate_norm, num_comps, num_trials, ensamble_method);
%[dred_factors_shuff, ~] = f_dred_train2(firing_rate_shuff, num_comps, num_trials, ensamble_method);

%f_dred_plot_factors(dred_data_list3,trial_types_dred, test_data_ind);


%% extract ensambles

%ens = f_extract_ensamble_from_scores(dred_factors1);
%ens = f_extract_ensamble_from_scores(dred_factors_shuff);



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

end