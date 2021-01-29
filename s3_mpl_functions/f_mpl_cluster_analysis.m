function f_mpl_cluster_analysis(data, ops)


for n_flip1 = 1:numel(ops.flip_to_analyze)
    for n_cond = 1:numel(ops.regions_to_analyze)
        n_flip = ops.flip_to_analyze(n_flip1);
        cond_name = ops.regions_to_analyze{n_cond};
        cdata = data.(cond_name);
        ctx_cells_mmn = cat(1,cdata.resp_cells_mmn_onset{:}) + cat(1,cdata.resp_cells_mmn_offset{:});
        ctx_traces = cat(1,cdata.trial_ave_mmn{:});
        % extract data
        if n_flip == 1
            ctx_col = [1 2 3];
        elseif n_flip == 2
            ctx_col = [4 5 6];
        elseif n_flip == 3
            ctx_col = [1 2 3 4 5 6];
        end
        
        ctx_cells_mmn = logical(sum(ctx_cells_mmn(:,ctx_col),2));
        ctx_traces = ctx_traces(:,:,ctx_col);
        ctx_traces = ctx_traces(ctx_cells_mmn,:,:);
        
        if ops.normalize_before
            norm_ctx_traces = IF_normalize_data(ctx_traces, 2);
        else
            norm_ctx_traces = ctx_traces;
        end    

        % need to reduce data dimensionality
        % not sure if max is the best way ,but ok

        if ops.dim_red_method == 1
            ctx_features = squeeze(max(norm_ctx_traces, [], 2));

%             % plot data
%             figure;
%             plot3(ctx_features(:,1), ctx_features(:,2), ctx_features(:,3), '.')
%             xlabel('control');
%             ylabel('redundant');
%             zlabel('deviant');
%             
%             for n_cell = 1:50
%                 figure;
%                 hold on;
%                 plot(ctx_traces(n_cell,:,1), 'k');
%                 plot(ctx_traces(n_cell,:,2), 'b');
%                 plot(ctx_traces(n_cell,:,3), 'r');
%                 title(n_cell);
%             end

        elseif ops.dim_red_method == 2
            % pca method
            disp('Dimensionality reduction with PCA...');

            % decompose each context separately
            p_coeff_c = cell(3,1);
            p_score_c = cell(3,1);
            p_explained_c = cell(3,1);
            p_mu_c = cell(3,1);

            for n_cntxt = 1:3          
                traces_ctx = norm_ctx_traces(:,:,n_cntxt);
                [p_coeff_c{n_cntxt}, p_score_c{n_cntxt} ,~,~,p_explained_c{n_cntxt} ,p_mu_c{n_cntxt}] = pca(traces_ctx');
            end


            if ops.extra_clustering_plots == 1
                % PCA info
                figure;
                hold on;
                for n_cntxt = 1:3
                    plot(p_explained_c{n_cntxt},ops.context_colors{n_cntxt});
                end
                title('Components variance explained');
                xlabel('Component');
                ylabel('Variance explained');
                legend(ops.context_name);
            end


            if ops.extra_clustering_plots == 1
                % reconstruct traces
                figure;
                hold on;
                for n_cntxt = 1:3
                    for n_comp = 1:ops.num_scores_to_use(n_cntxt)
                        plot(ops.analysis_t, p_score_c{n_cntxt}(:,n_comp), ops.context_colors{n_cntxt});
                    end
                end
                title('compontents');
                xlabel('time, sec')
                temp_legend = cell(sum(ops.num_scores_to_use),1);
                temp_index = 0;
                for ii = 1:numel(ops.num_scores_to_use)
                    for jj = 1:ops.num_scores_to_use(ii)
                        temp_index = temp_index + 1;
                        temp_legend{temp_index} = sprintf('%s, component %d, %.1f%%', ops.context_name{ii}, jj, p_explained_c{ii}(jj));
                    end
                end
                legend(temp_legend);
            end

%             % plot reconstructed
%             traces_rec = p_coeff(:,1:num_scores_use)*p_score(:,1:num_scores_use)' + p_mu';
%             figure;
%             hold on;
%             plot(ops.time_stim_window, traces_rec(3*plot_cell-2,:), 'k')
%             plot(ops.time_stim_window, traces_rec(3*plot_cell-1,:), 'b')
%             plot(ops.time_stim_window, traces_rec(3*plot_cell,:), 'r')
%             title(sprintf('After PCA trial average, cell %d, %d components', plot_cell, num_scores_use));
%             xlabel('time, sec');
%             legend('Control', 'Redundant', 'Deviant');


            if ops.extra_clustering_plots == 1
                for n_comp = 1:min(min(ops.num_scores_to_use),3)
                    figure;
                    plot3(p_coeff_c{1}(:,n_comp), p_coeff_c{2}(:,n_comp), p_coeff_c{3}(:,n_comp), '.k');
                    xlabel(ops.context_name{1});
                    ylabel(ops.context_name{2});
                    zlabel(ops.context_name{3});
                    title(sprintf('Principal component %d', n_comp));
                end
            end


            % create data with reduced dimenstions
            % (cell1[cont(comp1, comp2,...),red(...),dev(...)];cell2[...])
            ctx_features = [p_coeff_c{1}(:,1:ops.num_scores_to_use(1)).*rms(p_score_c{1}(:,1:ops.num_scores_to_use(1)),1),...
                                  p_coeff_c{2}(:,1:num_scores_to_use(2)).*rms(p_score_c{2}(:,1:ops.num_scores_to_use(2)),1),...
                                  p_coeff_c{3}(:,1:num_scores_to_use(3)).*rms(p_score_c{3}(:,1:ops.num_scores_to_use(3)),1)];

            ctx_traces_dim_red = zeros(size(norm_ctx_traces));
            for n_cntxt = 1:3                  
                ctx_traces_dim_red(:,:,n_cntxt) = p_coeff_c{n_cntxt}(:,1:ops.num_scores_to_use(1))*p_score_c{n_cntxt}(:,1:ops.num_scores_to_use(1))'+p_mu_c{n_cntxt}(:);
            end

            ctx_dim_red_error = squeeze(sum((norm_ctx_traces - ctx_traces_dim_red).^2,2));

            if ops.extra_clustering_plots == 1


                figure;
                plot(sum(ctx_dim_red_error,2));
                title('dimensionality reduction error LSE');


                % make a plot of component variances
                % first make legend
                temp_legend = cell(sum(ops.num_scores_to_use),1);
                if ops.num_scores_to_use(1) > 0
                    for ii = 1:ops.num_scores_to_use(1)
                        temp_legend{ii} = sprintf('cont %d', ii);
                    end
                end
                if ops.num_scores_to_use(2) > 0
                    for ii = (1:ops.num_scores_to_use(2)) + ops.num_scores_to_use(1)
                        temp_legend{ii} = sprintf('red %d', ii - ops.num_scores_to_use(1));
                    end
                end
                if ops.num_scores_to_use(3) > 0
                    for ii = (1:ops.num_scores_to_use(3)) + sum(ops.num_scores_to_use(1:2))
                        temp_legend{ii} = sprintf('dev %d', ii - sum(ops.num_scores_to_use(1:2)));
                    end
                end
                temp_legend = categorical(temp_legend);

                % variance within components
                figure;
                bar(temp_legend,var(ctx_features));
                title('Variance of rms power of each component');

                plot_cell = 1:2;
                for n_cell = plot_cell
                    % before dim reduction
        %             figure;
        %             plot(squeeze(ctx_traces(1,:,:)))

                    figure;
                    subplot(2,1,1);
                    hold on;
                    plot(ops.analysis_t, norm_ctx_traces(n_cell,:,1), 'k')
                    plot(ops.analysis_t, norm_ctx_traces(n_cell,:,2), 'b')
                    plot(ops.analysis_t, norm_ctx_traces(n_cell,:,3), 'r')
                    title(sprintf('Before PCA trial average, cell %d', n_cell));
                    xlabel('time, sec');
                    legend('Control', 'Redundant', 'Deviant');


                    % plot reconstructed again to make sure it works
                    subplot(2,1,2);
                    hold on;
                    plot(ops.analysis_t, ctx_traces_dim_red(n_cell,:,1), 'k')
                    plot(ops.analysis_t, ctx_traces_dim_red(n_cell,:,2), 'b')
                    plot(ops.analysis_t, ctx_traces_dim_red(n_cell,:,3), 'r')
%                     for n_cntxt = 1:3
%                         temp_trace = p_coeff_c{n_cntxt}(n_cell,1:params.num_scores_to_use(1))*p_score_c{n_cntxt}(:,1:params.num_scores_to_use(1))'+p_mu_c{n_cntxt}(n_cell);
% 
%                         plot(ops.analysis_t, temp_trace, ops.context_colors{n_cntxt});
%                     end
                    title(sprintf('After PCA, trial average, cell %d', n_cell));
                    xlabel('time, sec');

                    legend(sprintf('Control, %d comp, %.1f%%', ops.num_scores_to_use(1), sum(p_explained_c{1}(1:ops.num_scores_to_use(1)))),...
                           sprintf('Redundant %d comp, %.1f%%', ops.num_scores_to_use(1), sum(p_explained_c{2}(1:ops.num_scores_to_use(2)))),...
                           sprintf('Deviant %d comp, %.1f%%', ops.num_scores_to_use(1), sum(p_explained_c{3}(1:ops.num_scores_to_use(3)))));


                end
            end
        elseif ops.dim_red_method == 3
            % use means in larger time bins
            bin_size = 3;

%             % use decimate function
%             ctx_features = zeros(size(ctx_traces,1), round(size(ctx_traces,2)/bin_size), 3);
%             for ii = 1:size(ctx_traces,1)
%                 for jj = 1:3
%                     ctx_features(ii, :, jj) = decimate(ctx_traces(ii,:,jj),bin_size);
%                 end
%             end
%             figure;
%             hold on;
%             plot(1:33, ctx_traces(1, :, 1))
%             plot((1:3:33)+2, ctx_features(1, :, 1))
            % use mean
            ctx_features = zeros(size(norm_ctx_traces,1), round(size(norm_ctx_traces,2)/bin_size), 3);
            for ii = 1:floor(size(norm_ctx_traces,2)/bin_size)
                temp_window = (ii*bin_size-bin_size+1):(ii*bin_size);
                ctx_features(:,ii,:) = mean(norm_ctx_traces(:,temp_window,:),2);
            end
            clear temp_window;

            ctx_features = reshape(ctx_features, size(ctx_features,1), []);
%             figure;
%             hold on;
%             plot(1:33, ctx_traces(1, :, 1))
%             plot((1:3:33)+1, ctx_features(1, :, 1))

        end



%         % distribution, looks bit exponential, maybe take log
%         figure;
%         hold on
%         histogram((max_resp(:,1)))
%         histogram((max_resp(:,2)))
%         histogram((max_resp(:,3)))
%         title('distribution of responses');


        % have to normalize, according to jordan
        norm_ctx_features = IF_normalize_data(ctx_features, ops.norm_after_method);
%         if ops.norm_after_method == 1
%             % maybe no normalization 
%             norm_ctx_features = ctx_features;
%         elseif ops.norm_after_method == 2
%             norm_ctx_features = ctx_features./max(ctx_features, [], 2);
%         elseif ops.norm_after_method == 3
%             % center data around zero and normalize by their magnitude of vector
%             centered = ctx_features - mean(ctx_features,1);
%             norm_ctx_features = centered./sqrt(diag(centered*centered'));
%         elseif ops.norm_after_method == 4
%             % now normalize each of the axes by std
%             %max_resp_t = log(ctx_traces_dim_red);
%             max_resp_t = (ctx_features);
%             max_resp_t_std = std(max_resp_t);
%             norm_ctx_features = max_resp_t./max_resp_t_std;
%         end

        disp('Kmeans...');

%         % plot distributions of max responses
%         figure;
%         hold on;
%         for ii = 1:3
%             [f,xi] = ksdensity(norm_max_resp(:,ii));
%             plot(xi, f);
%         end
%         legend('Control', 'redundant', 'deviant');
%         title(sprintf('%s, %s distributions of normalized data', ops.conditions{n_cond}, ops.fig_title_run{n_flip}));
%         
%         plot_data = norm_max_resp;
%         
%         figure;
%         plot3(plot_data(:,1), plot_data(:,2), plot_data(:,3), '.k')
%         xlabel('Control');
%         ylabel('Redundant');
%         zlabel('Deviant');
%         title('Normalized scatter plot')
%         title(sprintf('%s, %s', ops.conditions{n_cond}, ops.fig_title_run{n_flip}));




        % cluster analysis
        k_repeats = 500;
        k_num_clusters = 8;



        [all_idx_data, all_C_data, all_sumd_data] = IF_run_k_means_clust(norm_ctx_features, k_num_clusters, k_repeats);

        data_dim = size(norm_ctx_features,2);

        % plot cluster centroids, tells how stable the clustering is
        C_best_num = cell(8,1);
        use_C_MLE = 1;

        for n_clust = ops.num_clusters(n_cond) %2:max(4,num_clusters(n_cond))

            all_centroids = all_C_data{n_clust};
            all_centroids_sorted = zeros(size(all_centroids));
            for ii = 1:size(all_centroids,3)
                [~, temp_ind] = sort(all_centroids(:,1,ii));
                all_centroids_sorted(:,:,ii) = all_centroids(temp_ind,:,ii);
            end
            all_centroids_flat = reshape(permute(all_centroids_sorted, [2 1 3]),[],k_repeats)';
            all_centroids_pooled = reshape(permute(all_centroids_sorted, [2 1 3]),data_dim,[])';

            % plot kernels of centroids
            if ops.extra_clustering_plots == 1
                figure;
                for ii = 1:min(5,data_dim)
                    subplot(min(5,data_dim),1,ii)
                    [f,xi] = ksdensity(squeeze(all_centroids_sorted(1,ii,:)), 'BandWidth', 0.001);
                    plot(xi,f);


                end
                suptitle('kde distribution of centroid x coordinate');
            end
            clear f xi all_centroids;

            % identify best cluster
            % here find the most common cluster


            some_dist_measure = sum(abs(all_centroids_flat).^4,2);
            C_MLE_num = find(some_dist_measure == mode(some_dist_measure));
            C_MLE = mode(all_centroids_sorted,3);

            if ops.extra_clustering_plots == 1
                if numel(unique(some_dist_measure)) > 1
                    cent_distr = hist(some_dist_measure, unique(some_dist_measure));
                else
                    cent_distr = 500;
                end
                figure;
                bar(cent_distr);
                title(sprintf('unique cluster centroids distribution, MLEdist = %.3f',sum(C_MLE(:).^4)));
            end
            clear cent_distr;


            % cluster the centroids and find closest one
            [~, C_of_C, ~] = kmeans(all_centroids_pooled, n_clust);

%             [~, temp_sort] = sort(C_of_C(:,1));
%             C_of_C_flat = reshape(C_of_C(temp_sort,:)',1,[]);
%             C_total_dist = sum((C_of_C_flat - all_centroids_flat).^4,2);

            C_nearest_num = find(some_dist_measure == mode(some_dist_measure));
            % in case there are identical solutions, pick the first
            C_nearest = all_centroids_sorted(:,:,C_nearest_num(1));

            if use_C_MLE
                C_best_num{n_clust} = C_MLE_num(1);
            else
                C_best_num{n_clust} = C_nearest_num(1);
            end
            clear C_nearest_num C_MLE_num all_centroids_sorted all_centroids_flat some_dist_measure;

            % some plots are sometimes good
            if ops.dim_red_method == 2
                if ops.centroid_plots == 1
                    for ii = 1:min(min(ops.num_scores_to_use),2)
                        temp_comp = [ii ops.num_scores_to_use(1)+ii ops.num_scores_to_use(1)+ops.num_scores_to_use(2)+ii];
                        figure;
                        hold on;
                        plot3(C_MLE(:,temp_comp(1)), C_MLE(:,temp_comp(2)), C_MLE(:,temp_comp(3)), 'or');
                        plot3(C_nearest(:,temp_comp(1)), C_nearest(:,temp_comp(2)), C_nearest(:,temp_comp(3)), 'og');
                        plot3(C_of_C(:,temp_comp(1)), C_of_C(:,temp_comp(2)), C_of_C(:,temp_comp(3)), '.r');
                        plot3(all_centroids_pooled(:,temp_comp(1)), all_centroids_pooled(:,temp_comp(2)), all_centroids_pooled(:,temp_comp(3)), '.k');
                        legend('Most likely centrold', 'Nearest data centroid to CofC', 'Centroids of centroids (CofC)', 'Centroids of data (CofD)')
                        title(sprintf('%s, %s %d clust, PC%d',cond_name, ops.fig_title_run{n_flip}, n_clust, ii));
                        xlabel(ops.context_name{1});
                        ylabel(ops.context_name{2});
                        zlabel(ops.context_name{3});
                    end
                end
            end
        end
        clear C_MLE C_nearest C_of_C all_centroids_pooled num_scores_to_use all_C_data;


        %C_final = all_C_data{params.num_clusters(n_cond)}(:,:,C_best_num{params.num_clusters(n_cond)}(1));
        icx_temp = all_idx_data(:,ops.num_clusters(n_cond),C_best_num{ops.num_clusters(n_cond)}(1));
        clear C_best_num all_idx_data;

        % reassing clusters so cluster sizes are decreasing with number
        icx_final = zeros(size(icx_temp));

        num_clust_cells = zeros(ops.num_clusters(n_cond),1);
        for n_clust = 1:ops.num_clusters(n_cond)
            num_clust_cells(n_clust) = sum(icx_temp == n_clust);
        end
        [~, icx_index_desc] = sort(num_clust_cells, 'desc');
        if sum(ops.custom_clust_order) > 0
            if numel(ops.custom_clust_order) == numel(icx_index_desc)
                icx_index_desc = icx_index_desc(ops.custom_clust_order);
            else
                warning('Adjust number of clusters in custom_clust_order');
            end
        end
        for n_clust = 1:ops.num_clusters(n_cond)
            icx_final(icx_temp == icx_index_desc(n_clust)) = n_clust;
        end
        clear icx_temp icx_index_desc num_clust_cells;



        % random data sets
        % generate shuffled datasets
        norm_ctx_features_shuff = zeros([size(norm_ctx_features),k_repeats]);

        for n_rep = 1:k_repeats
            for ii = 1:data_dim
                norm_ctx_features_shuff(:,ii, n_rep) = randsample(norm_ctx_features(:), size(norm_ctx_features,1));
            end
        end
%         figure;
%         bar(mean(mean(shuf_norm_max_resp,1),3));

        % generate norm distributed set
        norm_ctx_features_randn = randn([size(norm_ctx_features),k_repeats]);

        ctx_fet_std = std(norm_ctx_features);
        ctx_fet_mean = mean(norm_ctx_features);

        for ii = 1:data_dim
            norm_ctx_features_randn(:,ii,:) = ctx_fet_std(ii)*norm_ctx_features_randn(:,ii,:) + ctx_fet_mean(ii);
        end
        clear ctx_fet_std ctx_fet_mean data_dim;

        [~, ~, all_sumd_shuff] = IF_run_k_means_clust(norm_ctx_features_shuff, k_num_clusters, k_repeats);
        [~, ~, all_sumd_rand] = IF_run_k_means_clust(norm_ctx_features_randn, k_num_clusters, k_repeats);
        clear norm_ctx_features_randn norm_ctx_features_shuff norm_ctx_features;

        use_randn_shuff = 1;
        % normalize and compute mean and std
        all_sumdn_data = all_sumd_data/mean(all_sumd_data(1,:))*100;
        sumd_mean_data = mean(all_sumdn_data,2);
        %sumd_std_data = std(all_sumdn_data, [], 2);
        sumd_CI_data = prctile((all_sumdn_data - sumd_mean_data)', [95 5]);
        all_sumdn_shuff = all_sumd_shuff/mean(all_sumd_shuff(1,:))*100;
        sumd_mean_shuff = mean(all_sumdn_shuff,2);
        %sumd_std_shuff = std(all_sumdn_shuff, [], 2);
        sumd_CI_shuff = prctile((all_sumdn_shuff - sumd_mean_shuff)', [95 5]);
        all_sumdn_rand = all_sumd_rand/mean(all_sumd_rand(1,:))*100;
        sumd_mean_rand = mean(all_sumdn_rand,2);
        %sumd_std_rand = std(all_sumdn_rand, [], 2);
        sumd_CI_rand = prctile((all_sumdn_rand - sumd_mean_rand)', [95 5]);
        clear all_sumdn_data all_sumdn_shuff all_sumdn_rand all_sumd_shuff all_sumd_rand;


        %f plot stuff
        figure;
        hold on;
        errorbar(1:k_num_clusters, [0; diff(sumd_mean_data)],sumd_CI_data(1,:),sumd_CI_data(2,:), 'LineWidth', 3, 'color', [0 0 0]); %sumd_std_data/sqrt(k_repeats-1)
        if ~use_randn_shuff
            errorbar(1:k_num_clusters, [0; diff(sumd_mean_shuff)],sumd_CI_shuff(1,:),sumd_CI_shuff(2,:), 'LineWidth', 2);
        end
        errorbar(1:k_num_clusters, [0; diff(sumd_mean_rand)],sumd_CI_rand(1,:),sumd_CI_rand(2,:), 'LineWidth', 3, 'color', [0.5 0.5 0.5]);
        ylabel('change within-cluster distance');
        xlabel('num. clusters in solution');
        if use_randn_shuff
            legend('observed', 'shuffled');
        else
            legend('observed', 'shuffled', 'random');
        end
        title(sprintf('%s, %s', cond_name, ops.fig_title_run{n_flip}));
        temp_tick = get(gca, 'YTick');
        axis tight;
        ylim([min(temp_tick) max(temp_tick)]);
        clear temp_tick;
        IF_fig_ax_lim_minmax(gca);
        set(gca, 'YTick', [min(get(gca, 'YTick')), 0]);
        set(gca, 'XTick', [3 5 7]);
        clear sumd_mean_data sumd_mean_shuff sumd_mean_rand sumd_CI_data sumd_CI_shuff sumd_CI_rand;
        clear k_num_clusters k_repeats;


%         figure;
%         plot(sumd_means)
%         hold on
%         plot(sumd_shuf_means)
%         legend('Data', 'Shuffled');
%         title(sprintf('%s, %s', ops.conditions{n_cond}, ops.fig_title_run{n_flip}));



%         figure;
%         hold on;
%         for ii = 1:params.num_clusters(n_cond)
%             plot3(norm_ctx_traces(icx_final==ii,1), norm_ctx_traces(icx_final==ii,3), norm_ctx_traces(icx_final==ii,6), '.');
%         end
%         xlabel('Control');
%         ylabel('Redundant');
%         zlabel('Deviant');
%         title(sprintf('%s, %s: Clustering results', ops.conditions{n_cond}, ops.fig_title_run{n_flip}));



        % compute y_min and y_max for plots
        y_min = zeros(ops.num_clusters(n_cond),1);
        y_max = zeros(ops.num_clusters(n_cond),1);

        for n_clust = 1:ops.num_clusters(n_cond)
            traces_ave = mean(ctx_traces(icx_final==n_clust,:,:),1);
            traces_SEM = std(ctx_traces(icx_final==n_clust,:,:), [], 1)/sqrt(max(sum(icx_final==n_clust)-1,1));
            y_min(n_clust) = min(traces_ave(:) - traces_SEM(:));
            y_max(n_clust) = max(traces_ave(:) + traces_SEM(:));
        end
        y_min = min(y_min);
        y_max = max(y_max);

        % statistics, onset window average
        stat_all_onset_averages = squeeze(mean(ctx_traces(:,cdata.onset_window_frames{1},:),2));
        [~, stat_all_dev_cont_p, ~, stat_all_dev_cont] = ttest2(stat_all_onset_averages(:,3), stat_all_onset_averages(:,1));
        [~, stat_all_red_cont_p, ~, stat_all_red_cont] = ttest2(stat_all_onset_averages(:,2), stat_all_onset_averages(:,1));
        fprintf('-------Statistics-------\nPopulation average (%d cells), %s %s:\ndev - cont tstat = %.2f; p = %f\nred - cont tstat = %.2f; p = %f\n',size(ctx_traces,1), ops.conditions{n_cond}, ops.fig_title_run{n_flip}, stat_all_dev_cont.tstat,stat_all_dev_cont_p , stat_all_red_cont.tstat, stat_all_red_cont_p);
        clear stat_all_onset_averages;

        % plot average for all cells
        % 1 - control; 2 - red; 3 - dev
        traces_ave = squeeze(mean(ctx_traces,1));
        traces_SEM = squeeze(std(ctx_traces, [], 1))/sqrt(size(ctx_traces,1)-1);

        figure;
        hold on;
        for n_cntxt = 1:3
            % converting to z-scores is tricky, I dont do it here yet
            shadedErrorBar(cdata.trial_window{1}.trial_window_t, traces_ave(:,n_cntxt), traces_SEM(:,n_cntxt), 'lineprops', ops.context_colors{n_cntxt});
        end
        title(sprintf('Population average, %s %s', cond_name, ops.fig_title_run{n_flip}));  
        ylabel('z-score');
        xlabel(sprintf('time (sec)'));
        temp_lim = get(gca, 'YTick');
        axis tight;
        ylim([min([0, temp_lim]) max(temp_lim)]);
        clear temp_lim;
        IF_pop_figure_format(gca);


        disp('Cluster statistics')
        % plot averages of clusters
        for n_clust = 1:ops.num_clusters(n_cond)
            ctx_num_cells = sum(icx_final==n_clust);
            
            if ctx_num_cells>1
            % more statistics
                stat_clust_onset_averages = squeeze(mean(ctx_traces(icx_final==n_clust,cdata.onset_window_frames{1},:),2));
                [~, stat_clust_dev_cont_p, ~, stat_clust_dev_cont] = ttest2(stat_clust_onset_averages(:,3), stat_clust_onset_averages(:,1));
                [~, stat_clust_red_cont_p, ~, stat_clust_red_cont] = ttest2(stat_clust_onset_averages(:,2), stat_clust_onset_averages(:,1));
                fprintf('Cluster %d (%d cells), %s %s:\ndev - cont tstat = %.2f; p = %f\nred - cont tstat = %.2f; p = %f\n',n_clust ,ctx_num_cells, cond_name, ops.fig_title_run{n_flip}, stat_clust_dev_cont.tstat, stat_clust_dev_cont_p, stat_clust_red_cont.tstat, stat_clust_red_cont_p);
            else
                disp('No stat, 1 cell only')
            end

            traces_ave = squeeze(mean(ctx_traces(icx_final==n_clust,:,:),1));
            traces_SEM = squeeze(std(ctx_traces(icx_final==n_clust,:,:), [], 1))/sqrt(ctx_num_cells-1);
            %subplot(1,params.num_clusters(n_cond), n_clust);
            figure;

            % plot limited trials data
            %subplot(2,4,n_cond+(n_flip-1)*4)
            hold on;
            for n_cntxt = 1:3
                % converting to z-scores is tricky, I dont do it here yet
                shadedErrorBar(cdata.trial_window{1}.trial_window_t, traces_ave(:,n_cntxt), traces_SEM(:,n_cntxt), 'lineprops', ops.context_colors{n_cntxt});
            end
            axis tight;
            ylim([min([0,y_min]) y_max]);
            ylabel('z-score');
            xlabel('time (sec)');
            title(sprintf('cluster %d, %d cells; %s %s', n_clust, sum(icx_final==n_clust), cond_name, ops.fig_title_run{n_flip}));
            IF_pop_figure_format(gca);
        end
        clear y_min y_max traces_ave traces_SEM;


        [~, ctx_traces_sort_ind] = sort(icx_final);
        ctx_traces_sort = ctx_traces(ctx_traces_sort_ind,:,:);
        figure;
        for ii = 1:3
            subplot(1,3,ii);
            imagesc(cdata.trial_window{1}.trial_window_t,1:size(ctx_traces_sort,1),ctx_traces_sort(:,:,ii));
            colormap gray;
            caxis([0 6]);
            title(ops.context_name{ii});
            IF_fig_ax_lim_minmax(gca);
            if ii == 1
                xlabel('time (sec)');
                ylabel('cells');
            else
                set(gca, 'XTick', [], 'YTick', []);
            end
        end
        suptitle([cond_name ', ' ops.fig_title_run{n_flip} ': Sorted cells by clusters']);
        clear ctx_traces_sort ctx_traces_sort_ind;

%         % sort by time of first peak
%         max_time = zeros(size(ctx_traces,2),1);
%         for ii = 1:size(ctx_traces,2)
%             temp_max = find(ctx_traces(ii,:,3) == max(ctx_traces(ii,:,3)));
%             max_time(ii) = temp_max(1);
%         end
%         
%         [~, temp_ind] = sort(max_time);
%         temp_traces_sort = ctx_traces(temp_ind,:,:);
%         figure;
%         for ii = 1:3
%             subplot(1,3,ii);
%             imagesc(temp_traces_sort(:,:,ii))
%             colormap gray;
%             caxis([0 6]);
%             title(ops.context_name{ii});
%         end
%         suptitle([ops.conditions{n_cond} ', ' ops.fig_title_run{n_flip} ': Sorted cells by time of first peak']);
%                 if params.dim_red_method == 1


        % first get max or mean of traces
        ctx_max = max(ctx_traces, [], 2);
        %ctx_max = mean(ctx_traces, 2);

        clust_colors = {'r', 'g', 'k'};

        if ops.dim_red_method == 1
            % plot data
            figure;
            hold on;
            for n_clust = 1:ops.num_clusters(n_cond)
                plot3(ctx_max(icx_final==n_clust,1), ctx_max(icx_final==n_clust,2), ctx_max(icx_final==n_clust,3), '.', 'color', clust_colors{n_clust}, 'MarkerSize', 20);
            end
            xlabel('control');
            ylabel('redundant');
            zlabel('deviant');
            title('Max onset activity');
        end
        IF_fig_ax_lim_minmax(gca);
        clear ctx_traces ctx_max;

        % pie chart
        ctx_cells = zeros(ops.num_clusters(n_cond),1);
        for ii = 1:ops.num_clusters(n_cond)
            ctx_cells(ii) = sum(icx_final == ii);
        end
        figure;
        pie(ctx_cells);
        title('Totals');
        IF_figure_font_format(gca);

    end
end

end


%% internal functions

function [idx, C, sumd] = IF_run_k_means_clust(data, k_num_clusters, k_repeats)
    
    data_dim = size(data,2);

    idx = zeros(size(data,1), k_num_clusters, k_repeats);
    C = cell(k_num_clusters,1);
    sumd = zeros(k_num_clusters,k_repeats);
    % repeat many times

    if ndims(data) == 2
        data = repmat(data, [1 1 k_repeats]);
    end

    for n_clust = 1:k_num_clusters
        temp_C = zeros(n_clust,data_dim,k_repeats);
        temp_sumd = zeros(n_clust,k_repeats);
        for n_rep = 1:k_repeats
            [idx(:,n_clust,n_rep), temp_C(:,:,n_rep), temp_sumd(:,n_rep)] = kmeans(data(:,:,n_rep), n_clust);
        end

        C{n_clust} = temp_C;

        if size(temp_sumd,1) > 1
            temp_sumd = sum(temp_sumd,1);
        end
        sumd(n_clust,:) = temp_sumd;
    end
end

function norm_traces = IF_normalize_data(traces, norm_method)
    
    num_cells = size(traces,1);

    if norm_method == 1
        % maybe no normalization 
        norm_traces = traces;
    elseif norm_method == 2
        norm_traces = zeros(size(traces));
        for ii = 1:num_cells
            temp_cell_trace = traces(ii,:);
            norm_traces(ii,:) = traces(ii,:) / max(temp_cell_trace);
        end
    elseif norm_method == 3
        % center data around zero and normalize by their magnitude of vector
        centered = traces - mean(traces,1);
        % for each cell, normalize the feature vector magnitude
        norm_traces = centered./sqrt(diag(centered*centered'));
    elseif norm_method == 4
        % now normalize each cell by std of features
        %norm_traces = zeros(size(traces));
        %max_resp_t = log(ctx_traces_dim_red);
        max_resp_t = (traces);
        max_resp_t_std = std(max_resp_t, [], 2);
        norm_traces = max_resp_t./max_resp_t_std;
    end

%     figure;
%     plot(squeeze(norm_traces(3,:)))
%     hold on
%     plot(squeeze(traces(3,:)))
    
end

function IF_fig_ax_lim_minmax(fig_gca)

temp_tick = fig_gca.XTick;
fig_gca.XTick = [min(temp_tick) max(temp_tick)];
fig_gca.XMinorTick = 'on';
temp_tick = fig_gca.YTick;
fig_gca.YTick = [min(temp_tick), max(temp_tick)];
fig_gca.YMinorTick = 'on';
temp_tick = fig_gca.ZTick;
if sum(temp_tick) > 0
    fig_gca.ZTick = [min(temp_tick), max(temp_tick)];
    fig_gca.ZMinorTick = 'on';
end
    
IF_figure_font_format(fig_gca);
end

function IF_pop_figure_format(fig_gca)

fig_gca.XTick = [-0.1 0.6];
fig_gca.XMinorTick = 'on';
fig_gca.YTick = [0, max(fig_gca.YTick)];
fig_gca.YMinorTick = 'on';
temp_tick = fig_gca.ZTick;
if sum(temp_tick) > 0
    fig_gca.ZTick = [min(temp_tick), max(temp_tick)];
    fig_gca.ZMinorTick = 'on';
end

IF_figure_font_format(fig_gca);
end

function IF_figure_font_format(fig_gca)
fig_gca.FontSize = 10;
%fig_gca.FontWeight = 'bold';
fig_gca.FontWeight = 'normal';
fig_gca.TitleFontSizeMultiplier = 1;
fig_gca.LabelFontSizeMultiplier = 1;
fig_gca.TitleFontWeight = 'normal';
%fig_gca.LineWidth = 1;
end