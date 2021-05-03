function f_mpl_population_analysis_trials2(data, ops)
cv_data = cell(numel(ops.dred_params.trial_types_to_dred,1));
dim_est_st = struct('cond_name', [], 'n_dset', [], 'num_cells', [],...
    'num_cells_samp', [], 'num_comp_est', [], 'n_rep', []);
dd_idx = 1;
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    
    for n_dset = 1:numel(cdata.area)
        disp([cond_name, ' dset ' num2str(n_dset)]);
        
        trial_data_sort = cdata.trial_data_sort_wctx{n_dset};
        %trial_data_sort_sm = cdata.trial_data_sort_sm_wctx{n_dset};
        trial_types = cdata.trial_types_wctx{n_dset};
        trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
        
        
        %% cycle through analysis groups 
        for n_tt = 1:numel(ops.dred_params.trial_types_to_dred)
            %% select trials
            [tn_to_dred, trial_type_tag] = f_select_trial_type(ops.dred_params.trial_types_to_dred{n_tt}, cdata, n_dset, ops);
            tt_to_dred = ops.context_types_all(tn_to_dred);

            trials_idx_dred = logical(sum(trial_types == tt_to_dred(:)' ,2));
            trial_data_dred = trial_data_sort(:,:,trials_idx_dred);
            %trial_data_sm_dred = trial_data_sort_sm(:,:,trials_idx_dred);
            trial_peaks_dred = trial_peaks(:,trials_idx_dred);
            trial_types_dred = trial_types(trials_idx_dred);
            
            
            %% select responsive cells
            if ops.dred_params.use_responsive_cells
                resp_cells = cdata.peak_tuned_trials_combined{n_dset};
                resp_cells = logical(sum(resp_cells(:,tn_to_dred),2));
                trial_data_dred = trial_data_dred(resp_cells,:,:);
                %trial_data_sm_dred = trial_data_sm_dred(resp_cells,:,:);
                trial_peaks_dred = trial_peaks_dred(resp_cells,:);
            end
            [num_cells, ~, num_trials] = size(trial_data_dred);
            
            %%
            
            dr_params.cond_name = cond_name;
            dr_params.n_dset = n_dset;
            dr_params.trial_type_tag = trial_type_tag;
            dr_params.volume_period = cdata.proc_data{n_dset}.frame_data.volume_period;
            dr_params.tn_to_dred = tn_to_dred;
            dr_params.tt_to_dred = tt_to_dred;
            dr_params.trial_t = cdata.trial_window{n_dset}.trial_window_t;
            dr_params.ctx_mmn = ops.context_types_all(cdata.ctx_mmn{n_dset});
            dr_params.colors_clust = cat(2,ops.colors_list,ops.colors_list,ops.colors_list);
            
            
            %% hclustering
            num_clust = ops.dred_params.hclust.num_clust{n_tt};
            if ops.dred_params.do_hclust
                f_hcluster_trial(trial_peaks_dred, trial_types_dred, 'cosine', num_clust, dr_params, ops);
                f_hcluster_trial(trial_peaks_dred, trial_types_dred, 'ward', num_clust, dr_params, ops);
            end
            
            %% estimate dim with cross validation 
            if ops.dred_params.do_cv
                
                % randomize trial order but keep cell correlations
                rand_trial_indx = randperm(num_trials);     
                trial_data_trand = trial_data_dred(:,:,rand_trial_indx);
                %trial_types_dred_rand = trial_types_dred(rand_trial_indx);

                dr_params = f_make_dred_dir(dr_params, ops);

                dred_data_list = f_dim_red_cv(trial_data_trand, ops, dr_params);
                if isempty(cv_data{n_tt})
                    cv_data{n_tt} = struct();
                    cv_data{n_tt} = rmfield(dred_data_list,'dred_factors');
                else
                    cv_data{n_tt} = [cv_data{n_tt} rmfield(dred_data_list,'dred_factors')];
                end
            end
            
            %% estimate of dimensionality
            if ops.dred_params.do_dim_estimate
                trial_data_dred_sm = f_smooth_gauss(trial_data_dred, ops.ensemb.smooth_kernSD/cdata.proc_data{n_dset}.frame_data.volume_period);
                interval1 = 10;
                num_repeats = 10;
                dd_cells_range = [10:interval1:num_cells num_cells];
                for n_cellr = 1:numel(dd_cells_range)
                    for n_rep = 1:num_repeats
                        samp_idx = randsample(num_cells, dd_cells_range(n_cellr));
                        data_dim_est = f_ensemble_comp_data_dim(trial_data_dred_sm(samp_idx,:,:), 1);

                        %data_dim_est = f_ensemble_analysis_YS2(trial_data_sort_sm,trial_types_dred);
                        dim_est_st(dd_idx).cond_name = cond_name;
                        dim_est_st(dd_idx).n_dset = n_dset;
                        dim_est_st(dd_idx).tt_to_dred = tt_to_dred;
                        dim_est_st(dd_idx).trial_type_tag = trial_type_tag;
                        dim_est_st(dd_idx).num_cells = num_cells;
                        dim_est_st(dd_idx).num_cells_samp = dd_cells_range(n_cellr);
                        dim_est_st(dd_idx).num_comp_est = data_dim_est.num_comps;
                        dim_est_st(dd_idx).d_explained = data_dim_est.d_explained;
                        dim_est_st(dd_idx).n_rep = n_rep;
                        dim_est_st(dd_idx).var_thresh_prc = data_dim_est.var_thresh_prc;             
                        dd_idx = dd_idx+1;
                    end
                end
            end
            
            %% ensemble analysis
            if ops.dred_params.do_ensamble_analysis
                dr_params.trial_win_t = cdata.trial_window{n_dset}.trial_window_t;
                [~, dr_params.on_bin] = min(abs(ops.ensemb.onset_time-dr_params.trial_win_t));
                [~, dr_params.off_bin] = min(abs(ops.ensemb.offset_time-dr_params.trial_win_t));
                [~] = f_ensemble_analysis_YS2(trial_data_dred,trial_types_dred, dr_params, ops);
            end
            %if cdata.num_cells(n_dset) == max(cdata.num_cells) % strcmpi(cond_name, 'A2') %
                
            %end
            
        end
        
       
    end
end

% if ops.dred_params.do_cv
%     for n_tt = 1:numel(ops.dred_params.trial_types_to_dred)
%         cv_data2 = cv_data{n_tt}(strcmpi({cv_data{n_tt}.method}, 'svd'));
%         for n_cond = 1:numel(ops.regions_to_analyze)
%             cond_name = ops.regions_to_analyze{n_cond};
%             ccv_data = cv_data2(strcmpi({cv_data2.cond_name}, cond_name));
%             min_comp1 = zeros(max([ccv_data.n_dset]),1);
%             num_cells1 = zeros(max([ccv_data.n_dset]),1);
%             for n_dset = 1:max([ccv_data.n_dset])
%                 dccv_data = ccv_data([ccv_data.n_dset] == n_dset);
%                 num_fold = max([dccv_data.n_cv]);
%                 data1 = [[dccv_data.n_comp]', [dccv_data.train_err_sm]', [dccv_data.test_err_sm]'];
%                 data1 = reshape(data1,num_fold,[],3);
%                 mean_data = squeeze(mean(data1,1));
%                 [~, c_idx] = min(mean_data(:,3));
%                 min_comp1(n_dset) = mean_data(c_idx,1);
%                 num_cells1(n_dset) = dccv_data(1).num_cells;
%                 figure; hold on;
%                 plot(mean_data(:,1), mean_data(:,2))
%                 plot(mean_data(:,1), mean_data(:,3))
%             end
%             figure; plot(num_cells1, min_comp1, 'o')
% 
%         end
%     end
%     %save('cv_data_4_12_20', 'cv_data')
% end

if ops.dred_params.do_dim_estimate
    for n_tt = 1:numel(ops.dred_params.trial_types_to_dred)
        [~, trial_type_tag] = f_select_trial_type(ops.dred_params.trial_types_to_dred{n_tt}, cdata, n_dset, ops);
        dim_est_st2 = dim_est_st(strcmpi({dim_est_st.trial_type_tag}, trial_type_tag));
        colors1 = {[.5 .5 1], [1 .5 .5], [.5 1 .5]};
        colors2 = {[0 0 1], [1 0 0], [0 1 0]};
        figure; 
        s1 = subplot(1,2,1); hold on;
        s2 = subplot(1,2,2); hold on;
        pl1 = cell(numel(ops.regions_to_analyze),1);
        pl2 = cell(numel(ops.regions_to_analyze),1);
        dset_data_mean = cell(numel(ops.regions_to_analyze),1);
        for n_cond = 1:numel(ops.regions_to_analyze)
            cond_name = ops.regions_to_analyze(n_cond);
            dim_data1 = dim_est_st2(strcmpi({dim_est_st2.cond_name}, cond_name));
            for ii = 1:numel(dim_data1)
                dim_data1(ii).d_explained = sum(dim_data1(ii).d_explained);
            end
            dset_data_size = 0;
            for n_dset = 1:max([dim_data1.n_dset])
                dim_data2 = dim_data1([dim_data1.n_dset] == n_dset);
                data1 = [[dim_data2.num_cells_samp]' , [dim_data2.num_comp_est]', [dim_data2.d_explained]'];
                data2 = permute(mean(reshape(data1, max([dim_data2.n_rep]),[],3),1), [2 3 1]);
                dset_data_size = max([dset_data_size, size(data2,1)]);
            end
            dset_data = nan(dset_data_size-1,2,max([dim_data1.n_dset]));
            dset_data_lamb = nan(dset_data_size-1,2,max([dim_data1.n_dset]));
            for n_dset = 1:max([dim_data1.n_dset])
                dim_data2 = [dim_data1([dim_data1.n_dset] == n_dset)];
                data1 = [[dim_data2.num_cells_samp]' , [dim_data2.num_comp_est]', [dim_data2.d_explained]'];
                data2 = permute(mean(reshape(data1, max([dim_data2.n_rep]),[],3),1), [2 3 1]);
                subplot(s1)
                plot(data2(:,1),data2(:,2), 'Color', colors1{n_cond}, 'LineWidth', 0.1);
                subplot(s2)
                plot(data2(:,2),data2(:,3), 'Color', colors1{n_cond}, 'LineWidth', 0.1);
                dset_data(1:(size(data2,1)-1),:,n_dset) = data2(1:(size(data2,1)-1),1:2);
                dset_data_lamb(1:(size(data2,1)-1),:,n_dset) = data2(1:(size(data2,1)-1),[2, 3]);
            end
            dset_data_mean{n_cond} = nanmean(dset_data,3);
            dset_data_lamb_mean{n_cond} = nanmean(dset_data_lamb,3);
        end
        for n_cond = 1:numel(ops.regions_to_analyze)
            subplot(s1)
            pl1{n_cond} = plot(dset_data_mean{n_cond}(:,1),dset_data_mean{n_cond}(:,2), 'Color', colors2{n_cond}, 'LineWidth', 2);
            subplot(s2)
            pl2{n_cond} = plot(dset_data_lamb_mean{n_cond}(:,1),dset_data_lamb_mean{n_cond}(:,2), 'Color', colors2{n_cond}, 'LineWidth', 2);
        end
        subplot(s1)
        xlabel('Number of cells');
        ylabel('Estimated rank');
        %xlim([0 300]);
        %ylim([0 25]);
        legend([pl1{1} pl1{2} pl1{3}], ops.regions_to_analyze,'Location','southeast');
        title('Dimensionality of correlations', 'Interpreter', 'none');
        subplot(s2)
        xlabel('Data Rank');
        ylabel('Variance explained (%)');
        legend([pl2{1} pl2{2} pl2{3}], ops.regions_to_analyze,'Location','southeast');
        title('Variance in low rank data');
        suptitle(['trials analyzed ' trial_type_tag])
    end
end
end

