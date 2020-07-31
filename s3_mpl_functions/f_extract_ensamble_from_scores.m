function ens = f_extract_ensamble_from_scores(dred_factors, ens_thresh)

rectify_ens_score = 1;
plot_stuff = 1;
ens_min_size = 2;

method = dred_factors.method;

%num_cells = numel(dred_factors.dred_factors.means);

%%

[coeffs1, scores] = f_dred_get_coeffs(dred_factors);
[num_cells, num_comp] = size(coeffs1);
%%
%m_thresh = 0.2;
%z_thresh = 2;
ens = struct();
n_ens = 1;
for n_comp = 1:num_comp
    for n_thr = 1:numel(ens_thresh)
        if ens_thresh(n_thr)>0
            [coeffs_sort, sort_ind] = sort(coeffs1(:,n_comp), 'descend');
            cells = coeffs_sort > ens_thresh(n_thr);
            sign1 = 1;
        elseif ens_thresh(n_thr)<0
            [coeffs_sort, sort_ind] = sort(coeffs1(:,n_comp), 'ascend');
            cells = coeffs_sort < ens_thresh(n_thr);
            sign1 = -1;
        end
        if sum(cells) >= ens_min_size
            ens(n_ens).ens_cell_num = sort_ind(cells);
            ens(n_ens).cell_comp_mag = sign1*coeffs_sort(cells);
            ens(n_ens).n_comp = n_comp;
            ens(n_ens).ens_thresh = ens_thresh(n_thr);
            score_temp = sign1*scores(n_comp,:);
            if rectify_ens_score
                score_temp(score_temp<0) = 0;
            end
            ens(n_ens).ens_score = score_temp;
            ens(n_ens).num_cells = num_cells;
            ens(n_ens).method = method;
            n_ens = n_ens + 1;
        end
    end
end
%suptitle(method);


%%
if plot_stuff
    f1 = figure;
    f2 = figure;
    for n_comp = 1:num_comp
        [f, x] = ksdensity(coeffs1(:,n_comp), 'Bandwidth', 0.02);
        %[~, m_ind] = max(f);
        %k_mode = x(m_ind);
        %k_std = rms(coeffs1(:,n_comp));

        figure(f2); 
        subplot(ceil(num_comp/5),5,n_comp); hold on;
        plot(x,f);
        temp_ens = ens([ens.n_comp] == n_comp);
        for n_thr = 1:numel(temp_ens)
            thr_temp = temp_ens(n_thr).ens_thresh;
            line([thr_temp thr_temp], [0 max(f)], 'color', 'red', 'LineStyle', '--');
        end
        title(['comp ' num2str(n_comp)]);

        figure(f1);
        f_viz_subplot(ceil(num_comp/5),5,n_comp,1:num_cells,coeffs1(:,n_comp)); hold on;
        for n_thr = 1:numel(temp_ens)
            thr_temp = temp_ens(n_thr).ens_thresh;
            plot(ones(num_cells)*thr_temp, '--r')
        end
        plot(zeros(num_cells), '--g')
        title(num2str(n_comp));
    end
end

%             traceP = coeffsP'*trial_data_sm1_2d/(coeffsP'*coeffsP);  
%             cell1 = 2;
%             trace_cell_lr = traceP*ensP(n_comp).comp_mag(cell1);
%             trace_cell = trial_data_sm1_2d(ensP(n_comp).ens_cell_num(cell1),:);
%             figure; hold on;
%             plot(trace_cell)
%             plot(trace_cell_lr)
%             
%             trace_tca = reshape(dred_factors1.dred_factors.t_factors.U{2}(:,n_comp)*dred_factors1.dred_factors.t_factors.U{3}(:,n_comp)',1,[]);
%             figure; plot(traceP/max(traceP)); hold on;
%             plot(trace_tca/max(trace_tca))
% 
%             figure; plot(traceP/max(traceP)); hold on;
%             plot(trace_cell/max(trace_cell))

%%

% plot_cells = randsample(num_comps,example_plot);
% 
% if strcmpi(ensamble_method, 'PCA')
% 
%     %% PCA extracted components
%     % projecting columns(time bins) of Z(spiking raster) onto the
%             % population vector W (P = W*W')
%     % PCA extracted components
%     for n_pc = 1:num_comps
%         w = d_coeff(:,n_pc)/norm(d_coeff(:,n_pc)); % normalize
%         P = w*w';   % projection outer product  P = w*w'
%         P = P - diag(diag(P)); % set diagonal to zero so isolated activations do not contribute to R
% 
%         proj = P*firing_rate_norm; % projection of F onto w; proj = P*Zb;
%         
%         pop_vec(n_pc,:) = d_coeff(:,n_pc);
%         ens_act(n_pc,:) = d_score(:,n_pc);
%         
%         ens_cells{n_pc} = cell(2,1);
%         ens_cells{n_pc}{1} = find(d_coeff(:,n_pc)>(median(d_coeff(:,n_pc)) + 3*std(d_coeff(:,n_pc))));
%         ens_cells{n_pc}{2} = find(d_coeff(:,n_pc)<(median(d_coeff(:,n_pc)) - 3*std(d_coeff(:,n_pc))));
%         
%         
%         proj_dist = sum(proj);
%         Rb = diag(firing_rate_norm'*P*firing_rate_norm); % projection strength; Rb = Zb'*P*Zb (b is time bin)
%     end
%     
%     
% 
% elseif strcmpi(ensamble_method, 'AV')
%     %% assembly vector estimation method
% 
%     % project data to lower D space
%     corr_mat = firing_rate_norm*firing_rate_norm';
%     Pas = d_coeff(:,d_explained>data_thresh)*d_coeff(:,d_explained>data_thresh)';
%     N = Pas * corr_mat;        % project correlation matrix onto PC subspace
% 
%     % same as 
%     %FR2 = d_coeff(:,d_explained>data_thresh)*d_score(:,d_explained>data_thresh)';
%     %N2 = FR2*FR2';
% 
%     M = N'*N;                  % interaction matrix
%     
% %     corr_mat_shuff = firing_rate_cont_shuff*firing_rate_cont_shuff';
% %     Pas_shuff = s_coeff(:,d_explained>data_thresh)*s_coeff(:,d_explained>data_thresh)';
% %     N_shuff = Pas_shuff * corr_mat_shuff; 
% %     M_shuff = N_shuff'*N_shuff;
%     
%     % first need to binarize M
%     
%     
%     M_thresh = median(M(:))+3*std(M(:));
%     
%     %[idx,~,~]  = kmeans(M(:),2,'Replicates',5);
%     %M_thresh = max(min(M(idx==1)),min(M(idx==2)));
%     
%     M_binary = M>M_thresh;
%     
%     figure; imagesc(M)
%     figure; imagesc(M_binary)
% 
%     M_dist = 1-M_binary;
%     M_dist = M_dist - diag(diag(M_dist));
%     
%     Z = squareform(M_dist);
%     
%     Z = linkage(M_dist)
%     
%     figure;
%     [H,T,outperm] = dendrogram(Z);
%     
%     sum(T == 2)
%     
%     figure; imagesc(1-M_binary)
% 
% 
% elseif strcmp(ensamble_method, 'ICA')
%     %% ICA assembly estimation
%     n_pc = 1:num_comps;
%     w = d_coeff(:,n_pc)/norm(d_coeff(:,n_pc)); % normalize
%     P = w*w';   % projection outer product  P = w*w'
% 
%     P = P - diag(diag(P)); % set diagonal to zero so isolated activations do not contribute to R
% 
%     proj = P*firing_rate_norm;
% 
%     [icasig, A, W] = fastica(proj,'lastEig', num_comps); % 'lastEig' 
%     
%     for n_pc = 1:num_comps
%         pop_vec(n_pc,:) = A(:,n_pc);
%         ens_act(n_pc,:) = icasig(n_pc,:);
% 
%         ens_cells{n_pc} = cell(2,1);
%         ens_cells{n_pc}{1} = find(A(:,n_pc)>(median(A(:,n_pc)) + 3*std(A(:,n_pc))));
%         ens_cells{n_pc}{2} = find(A(:,n_pc)<(median(A(:,n_pc)) - 3*std(A(:,n_pc))));
%     end
% 
% 
% elseif strcmp(ensamble_method, 'NMF')
%     
%     plot_comp_num = 1:3;
%     
%     num_comp = 100;
%     % NMF minimizes D = norm(A-W*H,'fro')/sqrt(N*M)
%     % nuclear norm nucnorm = norm(svd(A),1);
%     [d_W,d_H] = nnmf(firing_rate,num_comp);
%     firing_rate_NMF = d_W*d_H;
%     
%     
%     comp_size = zeros(num_comp,1);
%     for n_comp = 1:num_comp
%         comp_size(n_comp) = sum(d_W(:,n_comp)>(max(d_W(:,n_comp))/2));
%     end
%     figure; plot(comp_size); title('comp sizes')
%     
%     fro_norm = zeros(num_comp,1);
%     nuc_norm = zeros(num_comp,1);
%     for n_comp = 1:num_comp
%         % compute norm
%         firing_rate_NMF1 = d_W(:,1:n_comp)*d_H(1:n_comp,:);
%         fro_norm(n_comp) = norm(firing_rate - firing_rate_NMF1, 'fro')/sqrt(num_cells*num_frames);
%         nuc_norm(n_comp) = norm(svd(firing_rate - firing_rate_NMF1),1)/sqrt(num_cells*num_frames);
%     end
%     %figure; plot(cumsum(fro_norm,'reverse'))
%     
%     figure; plot(fro_norm); title('Frobenius');
%     figure; plot(nuc_norm); title('Nuclear');
%     
%     
%     larger_clust = find(comp_size>1);
%     figure;
%     for n_cell_ind = 1:5
%         n_cell = larger_clust(n_cell_ind);
%         subplot(5,1,n_cell_ind); hold on;
%         plot(firing_rate(n_cell,:))
%         plot(firing_rate_NMF(n_cell,:))
%         title(['Cell ' num2str(n_cell) ' components = ' num2str(comp_size(n_cell))])
%     end
%     
%     
%     figure; imagesc(firing_rate); title('original')
%     figure; imagesc(firing_rate_NMF); title('NMF')
%     figure; imagesc(firing_rate - firing_rate_NMF); title('original - NMF')
%     
%     
%     
%     
%     
%     
%     int_comp = 89;
%     figure; plot(d_W(:,int_comp))
%     cells_in_comp = find(d_W(:,int_comp)>(max(d_W(:,int_comp))/2));
%     firing_rate_rank1 = d_W(:,int_comp)*d_H(int_comp,:);
%     for n_comp = int_comp
%         figure; imagesc(d_W(:,n_comp)*d_H(n_comp,:));
%         title(['Comp ' num2str(n_comp)]);
%     end
%     for n_comp = int_comp
%         figure; plot(d_W(:,n_comp))
%         title(['Comp ' num2str(n_comp)]);
%     end
%     for n_cell_ind = 1:numel(cells_in_comp)
%         n_cell = cells_in_comp(n_cell_ind);
%         figure; hold on;
%         plot(firing_rate(n_cell,:))
%         plot(firing_rate_rank1(n_cell,:));
%         title(['Cell ' num2str(n_cell)]);
%     end
%     
%     int_comp2 = find(comp_size>10);
%     figure; plot(d_W(:,int_comp2))
%     
% 
%     %figure; imagesc(d_H1*d_H2')
%     %figure; hold on; plot(d_H1(1,:));plot(d_H2(2,:))
%     
%     figure; imagesc(firing_rate); title('original')
%     figure; imagesc(firing_rate_NMF); title('NMF reduced')
%     figure; imagesc(firing_rate - firing_rate_NMF); title('remainder')
%     
%     norm(firing_rate, 'fro')/sqrt(num_cells*num_frames);
%     norm(firing_rate_NMF, 'fro')/sqrt(num_cells*num_frames);
%     norm(firing_rate - firing_rate_NMF, 'fro')/sqrt(num_cells*num_frames);
%     
%     n_cell = 1;
%     figure; hold on
%     plot(firing_rate(n_cell,:))
%     plot(firing_rate_NMF(n_cell,:))
%     title(['Cell ' num2str(n_cell)])
%     
%     [s_W,s_H] = nnmf(firing_rate_cont_shuffcirc,10);
%     
%     
%     
%     
%     
%     figure;
%     subplot(4,1,1); plot(d_H(plot_comp_num,:)');
%     title(['NNMF Components ' num2str(plot_comp_num)]);
%     subplot(4,1,2); plot(d_W(:,plot_comp_num));
%     subplot(4,1,3); plot(s_H(plot_comp_num,:)');
%     title(['NNMF shuffled Components ' num2str(plot_comp_num)]);
%     subplot(4,1,4); plot(s_W(:,plot_comp_num));
% 
% elseif strcmp(ensamble_method, 'SPCA')
%     %spca(X, Gram, K, delta, stop, maxSteps, convergenceCriterion, verbose)
%     % gram is X'*X covariance
%     % K is desired number of components
%     % delta is L2 ridge coefficient
%     % if delta is inf soft thresh is used (faster)
%     % stop is stopping criterion; 
%     %   stop positive is upper bound of L1 norm
%     %   stop negative is desired number of nonzero variables
%     %   stop 0 is regular PCA 
%     
%     K = 10;
%     delta = inf;
%     stop = 1000;
%     [SL, SD, L, D] = spca(firing_rate_norm', [], K, delta, stop);
%     
%     % SL is sparse loading vectors 
%     % SD is adjusted variances
%     % L and D loading variances of regular PCA
%     % paths loading path of each component as funct of iter number
%     
%     figure; plot(SL(:,5)); title('SL')
%     figure; plot(L(:,5)); title('L')
%     
%     
%     figure; plot(pop1(:,1))
%     
% end
% 




end