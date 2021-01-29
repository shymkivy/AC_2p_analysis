function f_mpl_cluster_analysis2(data, ops)

tn = [18 19 20; 28 29 30];
tt = ops.context_types_all(tn);

do_clust = 1;

pk_mag_all = cell(numel(ops.regions_to_analyze), numel(ops.flip_to_analyze));
pk_lat_all = cell(numel(ops.regions_to_analyze), numel(ops.flip_to_analyze));
for n_cond = 1:numel(ops.regions_to_analyze)
    for n_flip1 = 1:numel(ops.flip_to_analyze) 
        n_flip = ops.flip_to_analyze(n_flip1);
        cond_name = ops.regions_to_analyze{n_cond};
        cdata = data.(cond_name);
        ctx_cells_mmn = cat(1,cdata.peak_tuned_trials_combined_ctx{:});
        ctx_traces = cat(1,cdata.trial_ave_mmn{:});
        trial_window_t = cdata.trial_window{1}.trial_window_t;
        
        
        %% extract data
        if n_flip == 1
            ctx_col = [1 2 3];
            ctx_cells_mmn = logical(sum(ctx_cells_mmn(:,ctx_col),2));
            ctx_traces = ctx_traces(:,:,ctx_col);
            ctx_traces = ctx_traces(ctx_cells_mmn,:,:);
        elseif n_flip == 2
            ctx_col = [4 5 6];
            ctx_cells_mmn = logical(sum(ctx_cells_mmn(:,ctx_col),2));
            ctx_traces = ctx_traces(:,:,ctx_col);
            ctx_traces = ctx_traces(ctx_cells_mmn,:,:);
        elseif n_flip == 3
            ctx_col1 = [1 2 3];
            ctx_cells_mmn1 = logical(sum(ctx_cells_mmn(:,ctx_col1),2));
            ctx_traces1 = ctx_traces(:,:,ctx_col1);
            ctx_traces1 = ctx_traces1(ctx_cells_mmn1,:,:);
            ctx_col2 = [4 5 6];
            ctx_cells_mmn2 = logical(sum(ctx_cells_mmn(:,ctx_col2),2));
            ctx_traces2 = ctx_traces(:,:,ctx_col2);
            ctx_traces2 = ctx_traces2(ctx_cells_mmn2,:,:);
            ctx_traces = cat(1, ctx_traces1, ctx_traces2);
        end
        

        %%        
        [pk_mag, pk_latency] = max(ctx_traces,[],2);
        
        num_cells = size(pk_mag,1);
        
        pk_mag = squeeze(pk_mag);
        pk_mag_all{n_cond, n_flip} = pk_mag;
        pk_mag2 = (pk_mag);
        %pk_mag_norm = (pk_mag2 - mean(pk_mag2(:)))./std(pk_mag2(:));
        %pk_mag_norm = (pk_mag2 - mean(pk_mag2))./std(pk_mag2);
        pk_mag_norm = (pk_mag2)./std(pk_mag2);
        
        pk_latency = squeeze(pk_latency);
        pk_lat_all{n_cond, n_flip} = pk_latency;
        pk_latency_norm = (pk_latency - mean(pk_latency(:)))/std(pk_latency(:));
        
        figure;hold on;
        plot3(pk_mag_norm(:,1), pk_mag_norm(:,2), pk_mag_norm(:,3), '.');
        plot3(0,0,0, '.k')
        xlabel('cont'); ylabel('red'); zlabel('dev');
        title(sprintf('%s, flip%d', cond_name, n_flip));
        grid on;
        axis tight;
        ax1 = gca;
        lims1 = max([ax1.XLim; ax1.XLim; ax1.XLim]);
        xlim(lims1);
        ylim(lims1);
        zlim(lims1);
        
        
        %figure; plot(pk_mag(:,1), pk_mag(:,2), '.')
%         figure; hold on;
%         histogram(log(pk_mag(:,1)))
%         histogram(log(pk_mag(:,2)))
%         histogram(log(pk_mag(:,3)))
        
        if do_clust

            hc_params.num_clust = 1;
            hc_params.method = 'average';%'average', 'ward'
            hc_params.distance_metric = 'cosine';% 'cosine', 'euclidean, 'squaredeuclidean'
            hc_params.plot_dist_mat = 1;
            hc_params.plot_clusters = 0;
            hc_params.XY_label = 'Cells';
            hclust_out = f_hcluster_wrap(pk_mag_norm, hc_params);
            
            title(sprintf('%s flip%d, hclust', cond_name, n_flip));
            
            %% measure clust centers and distances to centers
            clusters1 = unique(hclust_out.clust_ident);
            num_clust = numel(clusters1);
            
            num_test_clust = 10;
            test_clust_mean_mean = zeros(num_test_clust,1);
            for n_clust = 1:num_test_clust
                clust_ident1 = cluster(hclust_out.Z, 'MaxClust', n_clust);
                clust_params = f_compute_clust_params(pk_mag_norm, clust_ident1);
                test_clust_mean_mean(n_clust) = mean(clust_params.clust_mean_dist);
            end
            
            
            %% now for shuff
            num_shuff_rep = 100;
            test_clust_mean_mean_shuff = zeros(num_test_clust, num_shuff_rep);
            for n_rep = 1:num_shuff_rep
                pk_mag_norm_shuff = zeros(size(pk_mag_norm));
                for ii = 1:3
                    pk_mag_norm_shuff(:,ii) = pk_mag_norm(randperm(num_cells),ii);
                end
                
                hc_params.plot_dist_mat = 0;
                hclust_out_shuff = f_hcluster_cell(pk_mag_norm_shuff, [], hc_params, ops);
                
                for n_clust = 1:num_test_clust
                    clust_ident1 = cluster(hclust_out_shuff.Z, 'MaxClust', n_clust);
                    clust_params = f_compute_clust_params(pk_mag_norm_shuff, clust_ident1);
                    test_clust_mean_mean_shuff(n_clust,n_rep) = mean(clust_params.clust_mean_dist);
                end
                
            end
            
            test_clust_mean_mean2 = test_clust_mean_mean/test_clust_mean_mean(1)*100;
            test_clust_mean_mean_shuff2 = test_clust_mean_mean_shuff/mean(test_clust_mean_mean_shuff(1,:))*100;
            
            figure; hold on;
            plot(test_clust_mean_mean2)
            errorbar(mean(test_clust_mean_mean_shuff2,2), std(test_clust_mean_mean_shuff2,[],2))
            
            figure; hold on;
            plot([0; diff(test_clust_mean_mean2)])
            errorbar([0; mean(diff(test_clust_mean_mean_shuff2),2)], [0; std(diff(test_clust_mean_mean_shuff2),[],2)])
            
            
            
            figure; hold on;
            plot(diff(test_clust_mean_mean))
            plot(diff(mean(test_clust_mean_mean_shuff,2)))
            
            test_clust_mean_mean
            
            
            Y = tsne(pk_mag_norm);
            %figure;
            %gscatter(Y(:,1),Y(:,2),hclust_out.clust_ident);
            plt_params = struct;
            plt_params.plot_mean = 0;
            plt_params.marker_type = '.';
            %plt_params.marker_size = 10;
            f_plot_comp_scatter(Y, hclust_out.clust_ident, plt_params);
            title(sprintf('%s flip%d, T-SNE', cond_name, n_flip));
            
            
            ax1 = cell(num_clust,1);
            ylim1 = [0 0];
            %ord1 = clust_ident(dend_order);
            for n_clust = 1:num_clust
                figure; hold on; axis tight;
                temp_traces = ctx_traces(hclust_out.clust_ident == n_clust,:,:);
                shadedErrorBar_YS(trial_window_t, mean(temp_traces(:,:,1),1), std(temp_traces(:,:,1),[],1)/max(sqrt(sum(hclust_out.clust_ident == n_clust)-1),1), ops.context_colors{1});
                shadedErrorBar_YS(trial_window_t, mean(temp_traces(:,:,2),1), std(temp_traces(:,:,2),[],1)/max(sqrt(sum(hclust_out.clust_ident == n_clust)-1),1), ops.context_colors{2});
                shadedErrorBar_YS(trial_window_t, mean(temp_traces(:,:,3),1), std(temp_traces(:,:,3),[],1)/max(sqrt(sum(hclust_out.clust_ident == n_clust)-1),1), ops.context_colors{3});
                ax1{n_clust} = gca;
                ylim1 = [min([ax1{n_clust}.YLim(1) ylim1(1)]) max([ax1{n_clust}.YLim(2) ylim1(2)])];
                title(sprintf('%s, flip%d; %.1f%s cells; clust %d',cond_name, n_flip, sum(hclust_out.clust_ident == n_clust)/numel(hclust_out.clust_ident)*100, '%',n_clust));
            end
            for n_ax = 1:numel(ax1)
                ax1{n_ax}.YLim = ylim1;
            end

            figure;hold on;
            for n_clust = 1:num_clust
                plot3(pk_mag_norm(hclust_out.clust_ident == n_clust,1), pk_mag_norm(hclust_out.clust_ident == n_clust,2), pk_mag_norm(hclust_out.clust_ident == n_clust,3), '.', 'Color', ops.colors_list{n_clust});
            end
            xlabel('cont'); ylabel('red'); zlabel('dev');
            grid on;
            axis tight;
            ax1 = gca;
            lims1 = max([ax1.XLim; ax1.XLim; ax1.XLim]);
            xlim(lims1);
            ylim(lims1);
            zlim(lims1);
            title(sprintf('%s, flip%d', cond_name, n_flip));
            %legend('1', '2', '3')
        
        end
        
        
%         pk_lat2 = pk_mag_norm.*pk_latency_norm;
%         pk_lat2 = (pk_lat2 - mean(pk_lat2(:)))/std(pk_lat2(:));
%         
%         
%         pk_mag_lat = [pk_mag_norm , pk_latency_norm];
%         pk_mag_lat2 = [pk_mag_norm , pk_lat2];
%         pk_mag_lat3 = [pk_mag_norm , pk_latency_norm(:,1:2)];
%         
%         SI_pk_mag = similarity_index(pk_mag, pk_mag);
%         SI_pk_mag_lat = similarity_index(pk_mag_lat, pk_mag_lat);
%         
%         figure;
%         imagesc(SI_pk_mag)
        

    end
end


% figure;hold on;
% for n_cond = 2:3
%     plot3(pk_mag_all{n_cond,3}(:,1), pk_mag_all{n_cond,3}(:,2), pk_mag_all{n_cond,3}(:,3), '.');
% end
% xlabel('cont'); ylabel('red'); zlabel('dev');
% title(sprintf('all cond flip%d', n_flip));
% legend('A1', 'A2', 'DF')
% grid on;
% axis tight;
% %xlim([0 20]);
% ylim([0 25]);
% zlim([0 25]);
% 
% figure; hold on
% histogram(pk_lat_all{1,3})
% histogram(pk_lat_all{2,3})
% histogram(pk_lat_all{3,3})
% 
% 
% figure;hold on;
% for n_flip = 1:2
%     plot3(pk_mag_all{3,n_flip}(:,1), pk_mag_all{3,n_flip}(:,2), pk_mag_all{3,n_flip}(:,3), '.');
% end
% xlabel('cont'); ylabel('red'); zlabel('dev');
% title(sprintf('DF'));
% legend('flip1', 'flip2')
% grid on;
% axis tight;
% %xlim([0 20]);
% ylim([0 25]);
% zlim([0 25]);
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