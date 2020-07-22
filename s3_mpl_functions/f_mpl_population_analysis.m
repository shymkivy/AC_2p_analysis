function f_mpl_population_analysis(data, ops) 
disp('Ensemble analysis...');
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    for n_dset = 1:cdata.num_dsets
        
        firing_rate_smooth = cat(1,cdata.firing_rate_smooth{n_dset,:});
        mmn_phase_mpl = cdata.proc_data{n_dset}.mmn_phase_mpl{1};
        firing_rate_cont = firing_rate_smooth(:,mmn_phase_mpl == 1);
        stim_frame_index = cdata.proc_data{n_dset}.stim_frame_index{1};
        

        %% my analysis
        % remove inactive cells
        
        active_cells = sum(firing_rate_cont,2) > 0;
        disp([num2str(sum(active_cells)) ' active cells']);
        
        firing_rate_cont(~active_cells,:) = [];
        
        out = f_ensemble_analysis_PCA_YS(firing_rate_cont);
        
        %%
        %%
        plot_method = 'tca';
        plot_num_comp = 10;
        dred_data_list2 = dred_data_list([dred_data_list.n_comp] == plot_num_comp);
        dred_data_list3 = dred_data_list2(strcmpi({dred_data_list2.method}, plot_method));
        
        %f_dred_plot_factors(dred_data_list3,trial_types_dred, test_data_ind);
        
        [~, trials_sort_idx] = sort(trial_types_dred, 'ascend');
        trial_data_sort_sort_act = trial_data_sort(active_cells,:,trials_sort_idx);
        trial_types_dred_sort = trial_types_dred(trials_sort_idx);
        trial_data_sm1 = f_smooth_gauss(trial_data_sort_sort_act, 200/volume_period);
        trial_data_sm1_2d = reshape(trial_data_sm1,num_cells_act,[]);
        [num_cells_act, num_frames] = size(trial_data_sort_sort_act);
        
        trial_data_sm1_tshuff = zeros(size(trial_data_sm1));
        for n_cell = 1:num_cells_act
            trial_data_sm1_tshuff(n_cell,:,:) = trial_data_sm1(n_cell,:,randsample(num_trials,num_trials));
        end
        trial_data_sm1_tshuff_2d = reshape(trial_data_sm1_tshuff,num_cells_act,[]);
        
        %figure; imagesc(trial_data_sm1_2d'*trial_data_sm1_2d)
        %figure; imagesc(trial_data_sm1_tshuff_2d'*trial_data_sm1_tshuff_2d)
        
%         trial_data_sm1_2d_shuff = zeros(size(trial_data_sm1_2d));
%         for n_cell = 1:num_cells_act
%             trial_data_sm1_2d_shuff(n_cell,:) = circshift(trial_data_sm1_2d(n_cell,:),ceil(rand(1)*num_frames));
%         end
%         trial_data_sm1_shuff = reshape(trial_data_sm1_2d_shuff,num_cells_act,[],num_trials);
%         
        [dred_factors1, ~] = f_dred_train(trial_data_sm1, plot_num_comp, plot_method);
        [dred_factors_shuff, ~] = f_dred_train(trial_data_sm1_tshuff, plot_num_comp, plot_method);
        
        
        
        %% get comp
        %tt_colors = [jet(10)];
        
        if ops.dred_params.dred_mmn
            c_map = [.8 .8 .8; .2 .6 1; 1 .2 .2;.6 .6 .6; .2 .6 1; 1 .2 .2];
        else
            c_map = jet(numel(tt_to_dred));
        end
        trial_types_dred_sort_colors = zeros(numel(trial_types_dred_sort),3);
        for n_tt = 1:numel(tt_to_dred)
            trial_types_dred_sort_colors(trial_types_dred_sort == tt_to_dred(n_tt),:) = trial_types_dred_sort_colors(trial_types_dred_sort == tt_to_dred(n_tt),:)+c_map(n_tt,:);
        end
        figure; imagesc(reshape(trial_types_dred_sort_colors,1,numel(trial_types_dred_sort),3));
        
        if strcmpi(plot_method, 'svd')
            coeffs1 = dred_factors1.dred_factors.coeffs;
        elseif strcmpi(plot_method, 'nmf')
            coeffs1 = dred_factors1.dred_factors.d_W;
        elseif strcmpi(plot_method, 'tca')
            coeffs1 = dred_factors1.dred_factors.t_factors.U{1};
            f_viz_ktensor_trials(dred_factors1.dred_factors.t_factors, trial_types_dred_sort)
            f_viz_ktensor_trials(dred_factors_shuff.dred_factors.t_factors, trial_types_dred_sort)
        elseif strcmpi(plot_method, 'fa')
            coeffs1 = dred_factors1.dred_factors.L;
        elseif strcmpi(plot_method, 'gpfa')
            coeffs1 = dred_factors1.dred_factors.estParams.C;
        end
        
        z_thresh = 2;
        figure;
        ensP = struct('ens_cell_num', [], 'comp_mag', [], 'num_cells', [], 'method', []);
        ensN = struct('ens_cell_num', [], 'comp_mag', [], 'num_cells', [], 'method', []);
        for n_comp = 1:plot_num_comp
            [f, x] = ksdensity(coeffs1(:,2));
            [~, m_ind] = max(f);
            k_mode = x(m_ind);
            k_std = rms(coeffs1(:,2)-k_mode);
            f_viz_subplot(2,5,n_comp,1:num_cells,coeffs1(:,n_comp))
            hold on;
            p_thresh = k_std*z_thresh+k_mode;
            n_thresh = -k_std*z_thresh+k_mode;
            plot(ones(num_cells)*p_thresh, '--r')
            plot(ones(num_cells)*n_thresh, '--r')
            plot(zeros(num_cells)+k_mode, '--g')
            title(num2str(n_comp));
            p_cells = coeffs1(:,n_comp) > p_thresh;
            ensP(n_comp).ens_cell_num = find(p_cells);
            ensP(n_comp).comp_mag = coeffs1(p_cells,n_comp);
            ensP(n_comp).num_cells = num_cells_act;
            ensP(n_comp).method = plot_method;
            %coeffsP = p_cells.*coeffs1(:,n_comp);
            
% 
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
%             
            
   
            n_cells = coeffs1(:,n_comp) < n_thresh;
            ensN(n_comp).ens_cell_num = find(n_cells);
            ensN(n_comp).comp_mag = coeffs1(n_cells,n_comp);
            ensN(n_comp).num_cells = num_cells_act;
            ensN(n_comp).method = plot_method;
        end
        suptitle(plot_method);
        
        
        for n_comp = 1:plot_num_comp
            
        end
        
        
        mean(coeffs1(:,1))
        
        %% NMF
        figure;
        for n_comp = 1:plot_num_comp
            subplot(plot_num_comp,1,n_comp);
            plot(dred_factors1.dred_factors.d_W(:,n_comp))
        end
        
        for n_comp = 1:plot_num_comp
            ens1 = dred_factors1.dred_factors.d_W(:,n_comp)>2*rms(dred_factors1.dred_factors.d_W(:,n_comp));
            f_plot_ensamble(ens1, trial_data_sort_act, trial_types_dred, n_comp)
        end
        
        
        %% tca
        t_fact1 = dred_factors1.dred_factors.t_factors;
        viz_ktensor(t_fact1);
        
        n_comp1 = 3;
        ens1 = t_fact1.U{1}(:,n_comp1)>(mean(t_fact1.U{1}(:,n_comp1))+3*std(t_fact1.U{1}(:,n_comp1)));
        ens2 = t_fact1.U{1}(:,n_comp1)<(mean(t_fact1.U{1}(:,n_comp1))-3*std(t_fact1.U{1}(:,n_comp1)));
       
        f_plot_ensamble(ens1, trial_data_sort_act, trial_types_dred, n_comp1)
        f_plot_ensamble(ens2, trial_data_sort_act, trial_types_dred, n_comp1)
        %%
      
        
        
        
        %% ensemble analysis jordan method
        
        %f_mpl_ensembleanalysis(firing_rate_smooth);
        
        
        %% ensemble stoixeion Luis method
%         disp('El Stoiceion...');
%         close all;
%         addpath('C:\Users\ys2605\Desktop\matlab\Stoixeion-master');
%         addpath('C:\Users\ys2605\Desktop\matlab\Stoixeion-master\AXIS');
%         addpath('C:\Users\ys2605\Desktop\matlab\Stoixeion-master\AGORA');
%         addpath('C:\Users\ys2605\Desktop\matlab\Stoixeion-master\PLEGADES');
%         addpath('C:\Users\ys2605\Desktop\matlab\Stoixeion-master\PLEGADES\DRAKONTOS');
%         
%         [Pools_coords] = Stoixeion(binary_firing_rate_smooth,[],traces_raw_cont);
        
        %% ensemble analysis with Luis method
        
        binary_firing_rate_smooth = double(firing_rate_cont>3*std(firing_rate_smooth,[],2));
        bin_no_act = sum(binary_firing_rate_smooth,2) == 0;
        binary_firing_rate_smooth(bin_no_act,:) = [];
        
        % first make minary data
        addpath('C:\Users\ys2605\Desktop\matlab\SVD_ensemble');
        params.pks = 4;
        params.ticut = 0.2;
        params.jcut = 0.06;
        [core_svd,state_pks_full,param] = findSVDensemble(binary_firing_rate_smooth,[],params);
        
        figure; hold on;
        temp_trace = core_svd{3};
        for n_cell = 1
            plot(sum(firing_rate_cont(temp_trace,:)))
        end
%         figure; hold on;
%         plot(firing_rate_smooth(10,:)*5)
%         plot(a(10,:))
        
        %% full trace first
        
       

        if ops.ensemb.select_upstates
            
        end
        

        if ops.ensemb.PCA_dim_reduction
            % center data for PCA
            c_mu = mean(firing_rate_cont,2);
            firing_rate_cont_cent = firing_rate_cont - c_mu;

            c_mu_shuff = mean(firing_rate_cont_shuffcirc,2);
            firing_rate_cont_cent_shuff = firing_rate_cont_shuffcirc - c_mu_shuff;

            %[U,S,V] = svd(firing_rate_freq_cen); % U == p_coeff

            % scores are eigentraces, so we have less data to use
            [d_coeff,d_score,~,~,d_explained,~] = pca(firing_rate_cont_cent');
            [s_coeff,s_score,~,~,s_explained,~] = pca(firing_rate_cont_cent_shuff');

            comp_to_use = sum(cumsum(d_explained)<ops.ensemb.pca_var_thresh);
            %comp_to_use = 50;

            figure; hold on;
            plot(d_explained);
            plot(s_explained);
            legend('var explained', 'shuffled');
            xlabel('component');
            ylabel('explained variance');
            title('PCA explained var');

            figure; hold on;
            plot(cumsum(d_explained))
            plot(cumsum(s_explained))
            plot(ones(num_cells,1)*ops.ensemb.pca_var_thresh, '--r')
            axis tight
            legend('var explained', 'shuffled', 'manual thresh');
            xlabel('cum components');
            ylabel('cumulative explained variance');
            title(sprintf('PCA cumulative explained var, using %d comp', comp_to_use));

            data_back = (d_score(:,1:comp_to_use)*d_coeff(:,1:comp_to_use)')'; % no mu necessary since was subtracted before

            figure; hold on;
            plot(firing_rate_cont_cent(1,:))
            plot(data_back(1,:))
            legend('original data', 'recreated from comp')

            trial_data_pc_sort = f_get_stim_trig_resp(d_score(:,1:50)', stim_frame_index(1:400), [1, 9]);
            trial_ave_pc = f_mpl_trial_average(trial_data_pc_sort,trial_types, 1:10, ops.baseline_removal_trial_ave);

            plot_pc = 1:10;
            figure;
            y_max = max(max(max(trial_ave_pc(plot_pc,:,:))));
            y_min = min(min(min(trial_ave_pc(plot_pc,:,:))));
            for n_trial = 1:10
                subplot(2,5,n_trial)
                plot(trial_ave_pc(plot_pc,:,n_trial)');
                ylim([y_min y_max]);
            end
        end
        
        
        
        
        %eva = evalclusters(x,clust,criterion)
        

    end
end

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


%% simulate some stupid data
        %sp_rast = randn(100, 10000);
        
        %[~,~,~,~,si_explained,~] = pca(sp_rast');
        %figure; imagesc(sp_rast*sp_rast'); title('RandomN spike raster SI');
        %figure; histogram(si_explained,10); title('RandomN spike raster PCA explained hist');

end


