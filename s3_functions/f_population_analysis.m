function f_population_analysis(data, ops) 
    for n_cond = 1:numel(ops.conditions_to_analyze)
        cond_name = ops.conditions{ops.conditions_to_analyze(n_cond)};

        context_resp_cells = data.(cond_name).context_resp_cells;
        trial_ave_z_resp = data.(cond_name).trial_ave_z(context_resp_cells(:,1),ops.win.no_base_window,:);

        context_mmn = data.(cond_name).context_mmn(context_resp_cells(:,1),:,:);
        num_resp_cells = size(trial_ave_z_resp,1);


        combine_traces = 1;
        normalize = 0;

        trial_ave_ctx = zeros(num_resp_cells, 3, ops.win.no_base_window_size);
        ctx_label = zeros(num_resp_cells,3);
        for n_cell = 1:num_resp_cells
            trial_ave_ctx(n_cell,:,:) = squeeze(trial_ave_z_resp(n_cell, :, context_mmn(n_cell,:,1)))';
            % normalization
            if normalize
                trial_ave_ctx(n_cell,:,:) = trial_ave_ctx(n_cell,:,:)./max(max(trial_ave_ctx(n_cell,:,:)));
            end
            ctx_label(n_cell, :) = 1:3; 
        end

        if combine_traces
            trial_ave_ctx_sort = (reshape(permute(trial_ave_ctx, [1 3 2]),size(trial_ave_ctx,1), []));
            %ctx_label_flat = 
        else
            trial_ave_ctx_sort = reshape(trial_ave_ctx, [], ops.win.no_base_window_size);
            ctx_label_flat = reshape(ctx_label, [], 1);


            % averages
            figure;
            hold on;
            plot(mean(trial_ave_ctx_sort(ctx_label_flat == 1,:)), 'k');
            plot(mean(trial_ave_ctx_sort(ctx_label_flat == 2,:)), 'b');
            plot(mean(trial_ave_ctx_sort(ctx_label_flat == 3,:)), 'r');
            title('Averages')
        end
        % flatten traces

        %pca
        [p_coeff_c, p_score_c ,~,~,p_explained_c ,p_mu_c] = pca(trial_ave_ctx_sort);
        num_scores = 1:3;
        rec_data = p_score_c(:,num_scores)*p_coeff_c(:,num_scores)' + p_mu_c;

        figure;
        plot(p_explained_c)
        title('percent explained')

        figure;
        hold on;
        plot(p_score_c(:,1))
        plot(p_score_c(:,2))
        plot(p_score_c(:,3))
        title('PCA scores pop-cov');

    %     figure;
    %     hold on;
    %     plot(mean(rec_data(ctx_label_flat == 1,:)), 'k');
    %     plot(mean(rec_data(ctx_label_flat == 2,:)), 'b');
    %     plot(mean(rec_data(ctx_label_flat == 3,:)), 'r');
    %     title('PCA rec pop-cov');


        [p_coeff_c2, p_score_c2 ,~,~,p_explained_c2 ,p_mu_c2] = pca(trial_ave_ctx_sort');

        num_scores = 1:3;
        rec_data2 = (p_score_c2(:,num_scores)*p_coeff_c2(:,num_scores)' + p_mu_c2)';

        figure;
        hold on;
        plot(p_score_c2(:,1))
        plot(p_score_c2(:,2))
        plot(p_score_c2(:,3))
        title('PCA scores time-cov');

    %     figure;
    %     hold on;
    %     plot(mean(rec_data2(ctx_label_flat == 1,:)), 'k');
    %     plot(mean(rec_data2(ctx_label_flat == 2,:)), 'b');
    %     plot(mean(rec_data2(ctx_label_flat == 3,:)), 'r');
    %     title('PCA rec time-cov');

        % do again but concatenate control-red-dev... maybe some pattern'



        figure;
        plot(trial_ave_ctx_sort(2,:))

        [p_coeff_c3, p_score_c3 ,~,~,p_explained_c3 ,p_mu_c3] = pca(trial_ave_ctx_sort);

        num_scores = 1:3;
        rec_data3 = p_score_c3(:,num_scores)*p_coeff_c3(:,num_scores)' + p_mu_c3;

        figure;
        hold on;
        plot(p_score_c3(:,1))
        plot(p_score_c3(:,2))
        plot(p_score_c3(:,3))
        title('PCA components pop-cov');


        figure;
        hold on;
        plot(p_coeff_c3(:,1))
        plot(p_coeff_c3(:,2))
        plot(p_coeff_c3(:,3))



        traject3 = reshape(p_coeff_c3(:,1:3),ops.win.no_base_window_size,[],3);


        figure;
        plot(traject3(:,1,3))

        figure;
        hold on;
        plot3(traject3(:,1,1),traject3(:,1,2),traject3(:,1,3), 'k', 'LineWidth', 2)
        plot3(traject3(:,2,1),traject3(:,2,2),traject3(:,2,3), 'b', 'LineWidth', 2)
        plot3(traject3(:,3,1),traject3(:,3,2),traject3(:,3,3), 'r', 'LineWidth', 2)
        legend('control', 'red', 'dev')
        xlabel('comp1');
        ylabel('comp2');
        zlabel('comp3');


        figure;
        plot(p_score_c3(:,2))

        figure;
        plot(p_score_c3(:,1), p_score_c3(:,2), 'o')

        figure;
        plot3(p_score_c3(:,1), p_score_c3(:,2),p_score_c3(:,3), 'o')

        score_to_analyze = 3;

        figure;
        plot(p_score_c3(:,score_to_analyze), 'o')

        [f,xi] = ksdensity(p_score_c3(:,score_to_analyze));
        figure;
        plot(xi, f)
        thresh = xi(find(f == max(f))) + 2*std(p_score_c3(:,score_to_analyze));
        clear f xi;

        pc1_cells = p_score_c3(:,score_to_analyze) > thresh;

        sum(pc1_cells);

        figure;
        hold on;
        plot(squeeze(mean(trial_ave_ctx(pc1_cells, 1, :))), 'k');
        plot(squeeze(mean(trial_ave_ctx(pc1_cells, 2, :))), 'b');
        plot(squeeze(mean(trial_ave_ctx(pc1_cells, 3, :))), 'r');



        figure;
        hold on;
        plot(mean(rec_data(ctx_label_flat == 1,:)), 'k');
        plot(mean(rec_data(ctx_label_flat == 2,:)), 'b');
        plot(mean(rec_data(ctx_label_flat == 3,:)), 'r');
        title('PCA rec pop-cov');

        [p_coeff_c4, p_score_c4 ,~,~,p_explained_c4 ,p_mu_c4] = pca(trial_ave_ctx_sort');
        num_scores = 1:3;
        rec_data4 = (p_score_c4(:,num_scores)*p_coeff_c4(:,num_scores)' + p_mu_c4)';

        figure;
        hold on;
        plot(p_score_c4(:,1))
        plot(p_score_c4(:,2))
        plot(p_score_c4(:,3))
        title('PCA components time-cov');

        figure;
        hold on;
        plot(mean(rec_data4), 'k');
        plot(mean(rec_data4(ctx_label_flat == 2,:)), 'b');
        plot(mean(rec_data4(ctx_label_flat == 3,:)), 'r');
        title('PCA rec time-cov');

    end

end