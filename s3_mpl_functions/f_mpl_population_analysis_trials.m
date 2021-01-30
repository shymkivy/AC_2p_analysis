function f_mpl_population_analysis_trials(data, ops)
cv_data = struct();
dim_est_st = struct('cond_name', [], 'n_dset', [], 'num_cells', [],...
    'num_cells_samp', [], 'num_comp_est', [], 'n_rep', []);
dd_idx = 1;
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    
    for n_dset = 1:numel(cdata.area)
        disp([cond_name, ' dset ' num2str(n_dset)]);
        
        stim_frame_index = cell(1,cdata.num_planes(n_dset));
        firing_rate = cell(1,cdata.num_planes(n_dset));
        trial_data = cell(1,cdata.num_planes(n_dset));
        for n_pl = 1:cdata.num_planes(n_dset)
            % plane specific params
            stim_frame_index{n_pl}= cdata.stim_frame_index{n_dset,n_pl};
            firing_rate{n_pl} = cdata.firing_rate{n_dset,n_pl};
            trial_data{n_pl} = f_get_stim_trig_resp(firing_rate{n_pl}, stim_frame_index{n_pl}, cdata.trial_num_baseline_resp_frames{n_dset});
        end
        firing_rate = cat(1,firing_rate{:});  
        trial_data = cat(1,trial_data{:});
        
        %% select trials
        trial_types = cdata.trial_types{n_dset};
        if ops.dred_params.dred_mmn
            if ops.dred_params.dred_mmn == 1
                tn_to_dred = cdata.ctx_mmn{n_dset}(1:3);
            elseif ops.dred_params.dred_mmn == 2
                tn_to_dred = cdata.ctx_mmn{n_dset}(4:6);
            elseif ops.dred_params.dred_mmn == 3
                tn_to_dred = cdata.ctx_mmn{n_dset};
            end
        else
            tn_to_dred = ops.dred_params.trial_types_to_dred;
        end
        %%
        
        %%
        tt_to_dred = ops.context_types_all(tn_to_dred);
        
        if strcmpi(ops.dred_params.trace_to_use, 'trials_specified')
            trials_idx_dred = logical(sum(trial_types == tt_to_dred' ,2));
            trial_data_dred = trial_data(:,:,trials_idx_dred);
            trial_types_dred = trial_types(trials_idx_dred);
        else
            trial_data_dred = trial_data;
            trial_types_dred = trial_types;
        end
%trial_data_2d = reshape(trial_data,num_cells,[]);
        


        %% select responsive cells
        if ops.dred_params.use_responsive_cells
            %resp_cells = logical(cdata.resp_cells_all_offset{n_dset}+cdata.resp_cells_all_onset{n_dset});
            resp_cells = cdata.peak_tuned_trials_combined{n_dset};
            resp_cells = logical(sum(resp_cells(:,tn_to_dred),2));
            trial_data_dred = trial_data_dred(resp_cells,:,:);
        end
        [num_cells, ~, num_trials] = size(trial_data_dred);
        
        
        %%
        if 1
            [trial_data_sort_wctx,trial_types_wctx] =  f_s3_add_ctx_trials(trial_data, trial_types, MMN_freq, ops);
            trials_idx_dred_wctx1 = logical(sum(trial_types_wctx == ops.context_types_all(20) ,2));
            trials_idx_dred_wctx2 = logical(sum(trial_types_wctx == ops.context_types_all(30) ,2));
            if (sum(trials_idx_dred_wctx1) > 15)
                f1 = figure;
                sp{1} = subplot(2,3,1);
                sp{2} = subplot(2,3,4);
                x = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag(resp_cells,trials_idx_dred_wctx1);
                [dend_order, clust_ident] = f_hierarch_clust(x', 5, sp);

                Y = tsne(x');
                subplot(2,3,3);
                gscatter(Y(:,1),Y(:,2),clust_ident);
                axis equal tight;
                title('T-SNE');
                
                figure(f1);
                sp{1} = subplot(2,3,2);
                sp{2} = subplot(2,3,5);
                x = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag(resp_cells,trials_idx_dred_wctx2);
                [dend_order, clust_ident] = f_hierarch_clust(x', 5, sp);

                Y = tsne(x');
                subplot(2,3,6);
                gscatter(Y(:,1),Y(:,2),clust_ident);
                axis equal tight;
                title('T-SNE');
                suptitle(sprintf('%s dset %d', cond_name, n_dset));
            end
        end
        
        %%
        %%
        % randomize trial order but keep cell correlations
        rand_trial_indx = randperm(num_trials);     
        trial_data_trand = trial_data_dred(:,:,rand_trial_indx);
        %trial_types_dred_rand = trial_types_dred(rand_trial_indx);
        
        
        %%
        dr_params.cond_name = cond_name;
        dr_params.n_dset = n_dset;
        dr_params.volume_period = cdata.proc_data{n_dset}.frame_data.volume_period;
        dr_params.tn_to_dred = tn_to_dred;
        dr_params.tt_to_dred = tt_to_dred;
        dr_params.trial_t = cdata.trial_window{n_dset}.trial_window_t;
        dr_params.ctx_mmn = ops.context_types_all(cdata.ctx_mmn{n_dset});
        if ops.dred_params.do_cv
            dred_data_list = f_dim_red_cv(trial_data_trand, ops, dr_params);
            if ~numel(fields(cv_data))
                cv_data = rmfield(dred_data_list,'dred_factors');
            else
                cv_data = [cv_data rmfield(dred_data_list,'dred_factors')];
            end
        end
        
        %% estimate of dimensionality
        
        trial_data_dred_sm = f_smooth_gauss(trial_data_dred, ops.ensemb.smooth_kernSD/cdata.proc_data{n_dset}.frame_data.volume_period);
        
        if ops.dred_params.dim_estimate
            interval1 = 10;
            num_repeats = 10;
            dd_cells_range = [10:interval1:num_cells num_cells];
            for n_cellr = 1:numel(dd_cells_range)
                for n_rep = 1:num_repeats
                    samp_idx = randsample(num_cells, dd_cells_range(n_cellr));
                    data_dim_est = f_ensemble_comp_data_dim(trial_data_dred_sm(samp_idx,:,:));

                    %data_dim_est = f_ensemble_analysis_YS2(trial_data_sort_sm,trial_types_dred);
                    dim_est_st(dd_idx).cond_name = cond_name;
                    dim_est_st(dd_idx).n_dset = n_dset;
                    dim_est_st(dd_idx).num_cells = num_cells;
                    dim_est_st(dd_idx).num_cells_samp = dd_cells_range(n_cellr);
                    dim_est_st(dd_idx).num_comp_est = data_dim_est.num_comps;
                    dim_est_st(dd_idx).d_explained = data_dim_est.d_explained;
                    dim_est_st(dd_idx).n_rep = n_rep;
                    dim_est_st(dd_idx).var_thresh_prc = data_dim_est.var_thresh_prc;
                    dim_est_st(dd_idx).trial_type_tag = ops.dred_params.trial_type_tag;
                    dd_idx = dd_idx+1;
                end
            end
        end
        
        %[~, trials_sort_idx] = sort(trial_types_dred, 'ascend');
        
        %% ensemble analysis
        
        dr_params.trial_win_t = cdata.trial_window{n_dset}.trial_window_t;
        [~, dr_params.on_bin] = min(abs(ops.ensemb.onset_time-dr_params.trial_win_t));
        [~, dr_params.off_bin] = min(abs(ops.ensemb.offset_time-dr_params.trial_win_t));
        
        %if cdata.num_cells(n_dset) == max(cdata.num_cells) % strcmpi(cond_name, 'A2') %
            %[~] = f_ensemble_analysis_YS2(trial_data_dred_sm,trial_types_dred, dr_params, ops);
        %end
        
        
        
        %%
%         plot_method = 'tca';
%         plot_num_comp = 10;
%         dred_data_list2 = dred_data_list([dred_data_list.n_comp] == plot_num_comp);
%         dred_data_list3 = dred_data_list2(strcmpi({dred_data_list2.method}, plot_method));
%         
%         %f_dred_plot_factors(dred_data_list3,trial_types_dred, test_data_ind);
%         
%         [~, trials_sort_idx] = sort(trial_types_dred, 'ascend');
%         trial_data_sort_sort_act = trial_data_sort(active_cells,:,trials_sort_idx);
%         trial_types_dred_sort = trial_types_dred(trials_sort_idx);
%         trial_data_sm1 = f_smooth_gauss(trial_data_sort_sort_act, 200/volume_period);
%         trial_data_sm1_2d = reshape(trial_data_sm1,num_cells_act,[]);
%         [num_cells_act, num_frames] = size(trial_data_sort_sort_act);
%         
%         trial_data_sm1_tshuff = zeros(size(trial_data_sm1));
%         for n_cell = 1:num_cells_act
%             trial_data_sm1_tshuff(n_cell,:,:) = trial_data_sm1(n_cell,:,randsample(num_trials,num_trials));
%         end
%         trial_data_sm1_tshuff_2d = reshape(trial_data_sm1_tshuff,num_cells_act,[]);
%         
        %figure; imagesc(trial_data_sm1_2d'*trial_data_sm1_2d)
        %figure; imagesc(trial_data_sm1_tshuff_2d'*trial_data_sm1_tshuff_2d)
        
%         trial_data_sm1_2d_shuff = zeros(size(trial_data_sm1_2d));
%         for n_cell = 1:num_cells_act
%             trial_data_sm1_2d_shuff(n_cell,:) = circshift(trial_data_sm1_2d(n_cell,:),ceil(rand(1)*num_frames));
%         end
%         trial_data_sm1_shuff = reshape(trial_data_sm1_2d_shuff,num_cells_act,[],num_trials);
%         
%         [dred_factors1, ~] = f_dred_train(trial_data_sm1, plot_num_comp, plot_method);
%         [dred_factors_shuff, ~] = f_dred_train(trial_data_sm1_tshuff, plot_num_comp, plot_method);
%         
%         
        
        %% get comp
        %tt_colors = [jet(10)];
        
%         if ops.dred_params.dred_mmn
%             c_map = [.8 .8 .8; .2 .6 1; 1 .2 .2;.6 .6 .6; .2 .6 1; 1 .2 .2];
%         else
%             c_map = jet(numel(tt_to_dred));
%         end
%         trial_types_dred_sort_colors = zeros(numel(trial_types_dred_sort),3);
%         for n_tt = 1:numel(tt_to_dred)
%             trial_types_dred_sort_colors(trial_types_dred_sort == tt_to_dred(n_tt),:) = trial_types_dred_sort_colors(trial_types_dred_sort == tt_to_dred(n_tt),:)+c_map(n_tt,:);
%         end
%         figure; imagesc(reshape(trial_types_dred_sort_colors,1,numel(trial_types_dred_sort),3));
% 
% 
%         
%         for n_comp = 1:plot_num_comp
%             
%         end
%         
%         
%         mean(coeffs1(:,1))
        
        %% NMF
%         figure;
%         for n_comp = 1:plot_num_comp
%             subplot(plot_num_comp,1,n_comp);
%             plot(dred_factors1.dred_factors.d_W(:,n_comp))
%         end
%         
%         for n_comp = 1:plot_num_comp
%             ens1 = dred_factors1.dred_factors.d_W(:,n_comp)>2*rms(dred_factors1.dred_factors.d_W(:,n_comp));
%             f_plot_ensamble(ens1, trial_data_sort_act, trial_types_dred, n_comp)
%         end
%         
        
        %% tca
%         t_fact1 = dred_factors1.dred_factors.t_factors;
%         viz_ktensor(t_fact1);
%         
%         n_comp1 = 3;
%         ens1 = t_fact1.U{1}(:,n_comp1)>(mean(t_fact1.U{1}(:,n_comp1))+3*std(t_fact1.U{1}(:,n_comp1)));
%         ens2 = t_fact1.U{1}(:,n_comp1)<(mean(t_fact1.U{1}(:,n_comp1))-3*std(t_fact1.U{1}(:,n_comp1)));
%        
%         f_plot_ensamble(ens1, trial_data_sort_act, trial_types_dred, n_comp1)
%         f_plot_ensamble(ens2, trial_data_sort_act, trial_types_dred, n_comp1)
%         
%         
        
        
%         dred_temp = dred_data_list(81:84);
%         fac1 = cell(4,1);
%         fac_sort1 = cell(4,1);
%         for n_cv = 1:4
%             fac1{n_cv} = dred_temp(n_cv).dred_factors.d_H;
%             fac_sort1{n_cv} = reshape(fac1{n_cv}, 1,[],sum(~test_data_ind(:,n_cv)));
%             trial_types2 = trial_types_dred(~test_data_ind(:,n_cv));
%             tt1 = unique(trial_types2);
%             figure;
%             for n_tt= 1:numel(tt1)
%                 subplot(2,3,n_tt)
%                 plot(squeeze(fac_sort1{n_cv}(:,:,trial_types2 == tt1(n_tt))), colors1{n_tt});
%                 title(num2str(tt1(n_tt)))
%             end
%             suptitle('component 1')
%         end
%         
%         
%         dred_temp = dred_data_list(177:180);
%         fac1 = cell(4,1);
%         fac_sort1 = cell(4,1);
%         figure; hold on;
%         for n_cv = 1:4
%             figure;
%             plot(dred_temp(n_cv).dred_factors.t_factors.U{2}) 
%         end
%         figure;
%         for n_cv = 1:4
%             trial_fac1 = dred_temp(n_cv).dred_factors.t_factors.U{3};
%             trial_types2 = trial_types_dred(~test_data_ind(:,n_cv));
%             subplot(2,2,n_cv); hold on;
%             for n_tr = 1:6
%                 plot(ones(numel(trial_fac1(trial_types2 == tt_to_dred(n_tr))),1)*n_tr, trial_fac1(trial_types2 == tt_to_dred(n_tr)), 'o', 'Color', colors1{n_tr})
%             end
%         end
%         
%         
%         for n_cv = 1:4   
%             fac1{n_cv} = dred_temp(n_cv).dred_factors.t_factors.U{2};
%             fac_sort1{n_cv} = reshape(fac1{n_cv}, 1,[],sum(~test_data_ind(:,n_cv)));
%             trial_types2 = trial_types_dred(~test_data_ind(:,n_cv));
%             tt1 = unique(trial_types2);
%             figure;
%             for n_tt= 1:numel(tt1)
%                 subplot(2,3,n_tt)
%                 plot(squeeze(fac_sort1{n_cv}(:,:,trial_types2 == tt1(n_tt))), colors1{n_tt});
%                 title(num2str(tt1(n_tt)))
%             end
%             suptitle('component 1')
%         end
%         
%         figure; plot(dred_data_list(1).dred_factors.coeffs)
%         
        
%         time_t = 1:n
%         
%         
%         
%         fac1 = reshape(dred_factors.dred_factors.scores
%         figure; plot(fac1)
%         figure; imagesc(squeeze(fac1)')
        
        
        
%         dred_params2 = ops.dred_params;
%         dred_params2.volume_period = volume_period;
%         dred_params2.save_path = [ops.file_dir ops.dred_params.save_path.(cond_name){n_dset}];
% 
%         f_dim_red_train(trial_data_sort(active_cells,:,:), dred_params2);
%         

%         temp_data = cat(1,data.(cond_name).firing_rate{n_dset,:});
%         num_cells = size(temp_data,1);
%         peaks = max(temp_data, [], 2);
%         
%         temp_data_norm = temp_data./peaks;
%         
%         peak_ind = zeros(num_cells,1);
%         for n_cell = 1:num_cells
%             temp_ind = find(temp_data_norm(n_cell,:) == 1);
%             peak_ind(n_cell) = temp_ind(1);
%         end
%         [~, sort_indx] = sort(peak_ind);
%         
%         %figure; imagesc(temp_data_norm)
%         figure;
%         ax1 = subplot(2,1,1);
%         imagesc(temp_data_norm(sort_indx,:))
%         title(sprintf('%s dataset %d sorted', cond_name, n_dset));
%         
%         ax2 = subplot(2,1,2); 
%         plot(mean(temp_data_norm,1));
%         axis tight;
%         linkaxes([ax1 ax2],'x');
    end
end

if ops.dred_params.do_cv
    
    cv_data = cv_data(strcmpi({cv_data.method}, 'svd'));
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze(n_cond);
        ccv_data = cv_data(strcmpi({cv_data.cond_name}, cond_name));
        min_comp1 = zeros(max([ccv_data.n_dset]),1);
        num_cells1 = zeros(max([ccv_data.n_dset]),1);
        for n_dset = 1:max([ccv_data.n_dset])
            dccv_data = ccv_data([ccv_data.n_dset] == n_dset);
            num_fold = max([dccv_data.n_cv]);
            data1 = [[dccv_data.n_comp]', [dccv_data.train_err_sm]', [dccv_data.test_err_sm]'];
            data1 = reshape(data1,num_fold,[],3);
            mean_data = squeeze(mean(data1,1));
            [~, c_idx] = min(mean_data(:,3));
            min_comp1(n_dset) = mean_data(c_idx,1);
            num_cells1(n_dset) = dccv_data(1).num_cells;
            %figure; hold on;
            %plot(mean_data(:,1), mean_data(:,2))
            %plot(mean_data(:,1), mean_data(:,3))
        end
        %figure; plot(num_cells1, min_comp1, 'o')
        
    end
    
    %save('cv_data_4_12_20', 'cv_data')
end

if ops.dred_params.dim_estimate
    colors1 = {[.5 .5 1], [1 .5 .5], [.5 1 .5]};
    colors2 = {[0 0 1], [1 0 0], [0 1 0]};
    f1 = figure; hold on;
    f2 = figure; hold on;
    pl1 = cell(numel(ops.regions_to_analyze),1);
    pl2 = cell(numel(ops.regions_to_analyze),1);
    dset_data_mean = cell(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze(n_cond);
        dim_data1 = dim_est_st(strcmpi({dim_est_st.cond_name}, cond_name));
        for ii = 1:numel(dim_data1)
            dim_data1(ii).d_explained = sum(dim_data1(ii).d_explained);
        end
        dset_data_size = 0;
        for n_dset = 1:max([dim_data1.n_dset])
            dim_data2 = dim_data1([dim_data1.n_dset] == n_dset);
            data1 = [[dim_data2.num_cells_samp]' , [dim_data2.num_comp_est]', [dim_data2.d_explained]'];
            data2 = squeeze(mean(reshape(data1, max([dim_data2.n_rep]),[],3),1));
            dset_data_size = max([dset_data_size, size(data2,1)]);
        end
        dset_data = nan(dset_data_size-1,2,max([dim_data1.n_dset]));
        dset_data_lamb = nan(dset_data_size-1,2,max([dim_data1.n_dset]));
        for n_dset = 1:max([dim_data1.n_dset])
            dim_data2 = [dim_data1([dim_data1.n_dset] == n_dset)];
            data1 = [[dim_data2.num_cells_samp]' , [dim_data2.num_comp_est]', [dim_data2.d_explained]'];
            data2 = squeeze(mean(reshape(data1, max([dim_data2.n_rep]),[],3),1));
            figure(f1)
            plot(data2(:,1),data2(:,2), 'Color', colors1{n_cond}, 'LineWidth', 0.1);
            figure(f2)
            plot(data2(:,2),data2(:,3), 'Color', colors1{n_cond}, 'LineWidth', 0.1);
            dset_data(1:(size(data2,1)-1),:,n_dset) = data2(1:(size(data2,1)-1),1:2);
            dset_data_lamb(1:(size(data2,1)-1),:,n_dset) = data2(1:(size(data2,1)-1),[2, 3]);
        end
        dset_data_mean{n_cond} = nanmean(dset_data,3);
        dset_data_lamb_mean{n_cond} = nanmean(dset_data_lamb,3);
    end
    for n_cond = 1:numel(ops.regions_to_analyze)
        figure(f1)
        pl1{n_cond} = plot(dset_data_mean{n_cond}(:,1),dset_data_mean{n_cond}(:,2), 'Color', colors2{n_cond}, 'LineWidth', 2);
        figure(f2)
        pl2{n_cond} = plot(dset_data_lamb_mean{n_cond}(:,1),dset_data_lamb_mean{n_cond}(:,2), 'Color', colors2{n_cond}, 'LineWidth', 2);
    end
    figure(f1)
    xlabel('Number of cells');
    ylabel('Estimated rank');
    xlim([0 300]);
    ylim([0 25]);
    legend([pl1{1} pl1{2} pl1{3}], ops.regions_to_analyze,'Location','southeast');
    title(['Dimensionality of correlations ' ops.dred_params.trial_type_tag], 'Interpreter', 'none');
    figure(f2)
    xlabel('Data Rank');
    ylabel('Variance explained (%)');
    legend([pl2{1} pl2{2} pl2{3}], ops.regions_to_analyze,'Location','southeast');
    title(['Variance in low rank data' ops.dred_params.trial_type_tag]);
end
end

