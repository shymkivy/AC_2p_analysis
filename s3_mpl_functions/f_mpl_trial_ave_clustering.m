function f_mpl_trial_ave_clustering(data, ops)

%tn = [18 19 20; 28 29 30];
tn = [3 4 5];
%tn = [20 30 18];
tt = ops.context_types_all(tn);

do_combined = 1;
do_normalize = 1;

num_tns = size(tn,1);

num_clust = 3;


pk_mag_all = cell(numel(ops.regions_to_analyze), num_tns+do_combined);
pk_lat_all = cell(numel(ops.regions_to_analyze), num_tns+do_combined);
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    trial_window_t = cdata.trial_window{1}.trial_window_t;
    
    if ops.use_zscores
        trial_ave = cat(1, cdata.trial_ave_z{:});
    else
        trial_ave = cat(1, cdata.trial_ave{:});
    end
    
    resp_cells = cat(1,cdata.peak_tuned_trials_combined{:});
    
    traces = cell(num_tns,1);
    tn_labels = cell(num_tns,1);
    for n_tn = 1:num_tns
        traces{n_tn} = trial_ave(logical(sum(resp_cells(:,tn(n_tn,:)),2)),:,tn(n_tn,:));
        tn_labels{n_tn} = num2str(tn(n_tn,:));
    end
    
    if do_combined
        traces1{1} = cat(1, traces{:});
        tn_labels1{1} =  ['combined ' num2str(tn(:)')];
    else
        traces1 = traces;
        tn_labels1 = tn_labels;
    end
    
    
    for n_tn = 1:numel(traces1)
        ctx_traces = traces1{n_tn};
        [~, ~, num_feat] = size(ctx_traces);
        plot_colors1 = ops.context_types_all_colors2(tn(n_tn,:));
        title_tag = sprintf('%s %s', cond_name, tn_labels1{n_tn});
        %%
        [pk_mag, pk_latency] = max(ctx_traces,[],2);
        pk_mag = squeeze(pk_mag);
        pk_latency = squeeze(pk_latency);
        pk_latency_t = trial_window_t(pk_latency);
        
        %% normalize mag
        if do_normalize
            %pk_mag_norm = (pk_mag2 - mean(pk_mag2(:)))./std(pk_mag2(:));
            %pk_mag_norm = (pk_mag - mean(pk_mag))./std(pk_mag);
            pk_mag_norm = (pk_mag)./std(pk_mag);
        else
            pk_mag_norm = pk_mag;
        end
        
        pk_mag_all{n_cond, n_tn} = pk_mag_norm;

        %%
        %pk_latency_norm = (pk_latency - mean(pk_latency(:)))/std(pk_latency(:));
        pk_lat_all{n_cond, n_tn} = pk_latency_t;
        
        %% plot peak latencies info
        latency_plot_index = logical((pk_latency_t > trial_window_t(1)).*(pk_latency_t < trial_window_t(end)));
        figure; 
        subplot(2,1,1);
        plot(pk_latency_t(latency_plot_index), pk_mag(latency_plot_index), '.');
        xlabel('time, sec');
        ylabel('resp mag');
        title([title_tag ' trial ave peak latencies'])
        subplot(2,1,2);
        histogram(pk_latency_t(latency_plot_index))
        xlabel('time, sec');
        ylabel('count');
        
        %% plot scatter before
        plt_params = struct;
        plt_params.plot_mean = 0;
        plt_params.marker_type = '.';
        plt_params.marker_size = 10;
        plt_params.equalize_lims = 1;
        f_plot_comp_scatter(pk_mag_norm, [], plt_params);
        title(title_tag);
        xlabel(ops.context_types_labels{tn(n_tn,1)});
        ylabel(ops.context_types_labels{tn(n_tn,2)});
        
        %% cluster
        hc_params.num_clust = num_clust;
        hc_params.estimate_clust_num = 1;
        hc_params.method = 'average'; %'average', 'ward', ''
        hc_params.distance_metric = 'cosine'; % 'cosine', 'euclidean, 'squaredeuclidean' hammilarity
        hc_params.plot_dist_mat = 1;
        hc_params.plot_clusters = 1;
        hc_params.XY_label = 'Cells';
        hc_params.title_tag = title_tag;
        hclust_out = f_hcluster_wrap(pk_mag_norm, hc_params);


        %% rearange the clust number for better coloring
        % first get mean vectors for each clust
        men_cl_vec = zeros(num_clust, num_feat);
        for n_cl = 1:num_clust
            men_cl_vec(n_cl,:) = mean(pk_mag_norm(hclust_out.clust_ident == n_cl,:));
        end
        
        align_mat = f_align_vec(diag(ones(num_feat,1)), men_cl_vec, hc_params.metric);
        
        clust_ident_mod = align_mat(hclust_out.clust_ident,1);
        
        %% plot scatte after
        f_plot_comp_scatter(pk_mag_norm, clust_ident_mod, plt_params, plot_colors1);
        title(title_tag);
        xlabel(ops.context_types_labels{tn(n_tn,1)});
        ylabel(ops.context_types_labels{tn(n_tn,2)});
        
        f_plot_comp_scatter(pk_mag, clust_ident_mod, plt_params, plot_colors1);
        title(title_tag);
        xlabel(ops.context_types_labels{tn(n_tn,1)});
        ylabel(ops.context_types_labels{tn(n_tn,2)});
        
        %% tsne
        Y = tsne(pk_mag_norm);
        %figure;
        %gscatter(Y(:,1),Y(:,2),hclust_out.clust_ident);
        plt_params.equalize_lims = 0;
        f_plot_comp_scatter(Y, clust_ident_mod, plt_params, plot_colors1);
        title(sprintf('%s, T-SNE',title_tag));

        %% trial cell averages for each clust
        ax1 = cell(num_clust,1);
        ylim1 = [0 0];
        %ord1 = clust_ident(dend_order);
        for n_clust = 1:num_clust
            figure; hold on; axis tight;
            clust_indx = clust_ident_mod == n_clust;
            temp_traces = ctx_traces(clust_ident_mod == n_clust,:,:);
            shadedErrorBar_YS(trial_window_t, mean(temp_traces(:,:,1),1), std(temp_traces(:,:,1),[],1)/max(sqrt(sum(clust_indx)-1),1), ops.context_colors{1});
            shadedErrorBar_YS(trial_window_t, mean(temp_traces(:,:,2),1), std(temp_traces(:,:,2),[],1)/max(sqrt(sum(clust_indx)-1),1), ops.context_colors{2});
            shadedErrorBar_YS(trial_window_t, mean(temp_traces(:,:,3),1), std(temp_traces(:,:,3),[],1)/max(sqrt(sum(clust_indx)-1),1), ops.context_colors{3});
            ax1{n_clust} = gca;
            ylim1 = [min([ax1{n_clust}.YLim(1) ylim1(1)]) max([ax1{n_clust}.YLim(2) ylim1(2)])];
            title(sprintf('%s; %.1f%s cells; clust %d',title_tag, sum(clust_indx)/numel(clust_ident_mod)*100, '%',n_clust));
        end
        for n_ax = 1:numel(ax1)
            ax1{n_ax}.YLim = ylim1;
        end

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
