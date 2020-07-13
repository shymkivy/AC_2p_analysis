function data = f_mpl_compute_tuning(data, ops)

disp('Computing tunning...');

% load_file = 0;
% stat_save = struct;
% if exist([ops.file_dir '\' ops.paradigm_type '_' ops.stat.save_est_samp '.mat'], 'file')
%     load([ops.file_dir '\' ops.paradigm_type '_' ops.stat.save_est_samp '.mat'], 'stat_save');
%     if numel(fields(stat_save))
%         load_file = 1;
%     end
% end
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
%     if ~load_file
%         %stat_save.(cond_name).onset_thresh = cell(data.(cond_name).num_dsets,1);
%         %stat_save.(cond_name).offset = cell(data.(cond_name).num_dsets,1);
%     end
    for n_dset = 1:data.(cond_name).num_dsets        
        cdata = data.(cond_name);
        
        trial_data_sort = cell(1,cdata.num_planes(n_dset));
        trial_num_baseline_resp_frames = cdata.trial_num_baseline_resp_frames{n_dset};
        for n_pl = 1:cdata.num_planes(n_dset)
            stim_frame_index = cdata.stim_frame_index{n_dset,n_pl};
            firing_rate = cdata.firing_rate_smooth{n_dset,n_pl};
            trial_data_sort{n_pl} = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_num_baseline_resp_frames); 
        end
        trial_data_sort = cat(1,trial_data_sort{:});
        trial_types = cdata.trial_types{n_dset};
        
        if ops.remove_early_dev
            for n_tr = 2:numel(trial_types)
                if trial_types(n_tr) == 170
                    if (trial_types(n_tr-1) < 103) || (trial_types(n_tr-1) == 170)
                        trial_types(n_tr) = 0;
                    end
                end
                if trial_types(n_tr) == 270
                    if (trial_types(n_tr-1) < 203) || (trial_types(n_tr-1) == 270)
                        trial_types(n_tr) = 0;
                    end
                end
            end
            data.(cond_name).trial_types{n_dset} = trial_types;
        end

        [trial_data_sort_pr,trial_types_pr] =  f_add_red_pool_trials(trial_data_sort, trial_types, ops);
        data.(cond_name).trial_data_sort_pr{n_dset,1} = trial_data_sort_pr;
        data.(cond_name).trial_types_pr{n_dset,1} = trial_types_pr;
        
        MMN_freq = cdata.MMN_freq{n_dset};
        dset_params.cond_name = cond_name;
        dset_params.n_dset = n_dset;
        dset_params.trial_window_t = cdata.trial_window_t{n_dset};
        dset_params.onset_window_frames = cdata.onset_window_frames{n_dset};
        dset_params.offset_window_frames = cdata.offset_window_frames{n_dset};
        dset_params.trial_types = trial_types_pr;
        dset_params.MMN_freq = MMN_freq;
        
%% create MMN sorted files
        
        % full set is flipped, no need to adjust this
        %ctx_mmn = [MMN_freq(2), 10 + ops.redundent_to_analyze, 19, ...
        %           MMN_freq(1), 19 + ops.redundent_to_analyze, 28];             
        ctx_mmn = [MMN_freq(2), 19, 20, ...
                   MMN_freq(1), 29, 30];  
        
        dset_params.ctx_mmn = ctx_mmn;

        %% get trial averages compute tunning
         sig_thresh = [];        
%         if load_file
%             if isfield(stat_save.(cond_name), 'sig_thresh_all')
%                 sig_thresh_all = stat_save.(cond_name).sig_thresh_all{n_dset};
%             end
%         end        

        tuning_all = f_get_tuning_dset(trial_data_sort_pr, ops.context_types_all, sig_thresh, dset_params, ops);
        data.(cond_name).tuning_all{n_dset,1} = tuning_all;

        %%
        z_thresh_all = max(mean(tuning_all.trace_tuning.stat_trace.z_factors,3),[],2);
        
        trial_ave = f_mpl_trial_average(trial_data_sort_pr,trial_types_pr, ops.context_types_all, 'none');
        trial_ave_z = trial_ave./z_thresh_all;
        %trial_ave = (trial_ave-means_all)./z_thresh_all;


        data.(cond_name).trial_ave{n_dset,1} = trial_ave;
        data.(cond_name).trial_ave_z{n_dset,1} = trial_ave_z;
        data.(cond_name).trial_ave_mmn{n_dset,1} = trial_ave(:,:,ctx_mmn);
        data.(cond_name).trial_ave_mmn_z{n_dset,1} = trial_ave_z(:,:,ctx_mmn);
        %data.(cond_name).resp_cells_mmn_onset{n_dset} = data.(cond_name).resp_cells_all_onset{n_dset}(:,ctx_mmn);
        %data.(cond_name).resp_cells_mmn_offset{n_dset} = data.(cond_name).resp_cells_all_offset{n_dset}(:,ctx_mmn);
        data.(cond_name).ctx_mmn{n_dset,1} = ctx_mmn;
        data.(cond_name).peak_tuned_trials_onset{n_dset,1} = tuning_all.peak_tuned_trials_onset;
        data.(cond_name).peak_tuned_trials_onset_ctx{n_dset,1} = tuning_all.peak_tuned_trials_onset_ctx;
        data.(cond_name).peak_tuned_trials_offset{n_dset,1} = tuning_all.peak_tuned_trials_offset;
        data.(cond_name).peak_tuned_trials_offset_ctx{n_dset,1} = tuning_all.peak_tuned_trials_offset_ctx;
        data.(cond_name).peak_tuned_trials_full{n_dset,1} = tuning_all.peak_tuned_trials_full;
        data.(cond_name).peak_tuned_trials_full_ctx{n_dset,1} = tuning_all.peak_tuned_trials_full_ctx;
        data.(cond_name).peak_tuned_trials_combined{n_dset,1} = tuning_all.peak_tuned_trials_combined;
        data.(cond_name).peak_tuned_trials_combined_ctx{n_dset,1} = tuning_all.peak_tuned_trials_combined_ctx;


%% deviance trial responses

%         trial_ave_freq_long = f_mpl_trial_average(trial_data_sort_long,trial_types, [170, 270], 'none');
%     
%         figure; hold on;
%         plot(trial_window_t_long,squeeze(mean(trial_ave_freq_long(:,:,1))))
%         plot(trial_window_t_long,squeeze(mean(trial_ave_freq_long(:,:,2))))
%         legend('MMN 170', 'flip MMN 270')
%         axis tight;
        
%         if ops.stat.plot_examples
%             for n_cell_ind = 1:numel(plot_cells)%900:910
%                 n_cell = plot_cells(n_cell_ind);
%                 
%                 figure;
%                 plot(trial_ave_freq_long(n_cell,:,1))
%                 title(sprintf('Cell %d', n_cell));
%             end
%         end
%     
%         num_cells = size(firing_rates,1);
%         resp_cells_zmag = data.(cond_name).resp_cells_zmag{n_dset};
%         resp_cells_zmag2 = resp_cells_zmag(:,1:10);
%         
%         resp_indx = zeros(num_cells,1);
%         for n_cell = 1:num_cells
%             temp_trace = resp_cells_zmag2(n_cell,:);
%             [~, temp_indx] = find(max(temp_trace) == temp_trace);
%             resp_indx(n_cell) = temp_indx(1);
%         end
%         
%         [~, resp_indx_sorted] = sort(resp_indx);
%         
%         figure; imagesc(resp_cells_zmag2(resp_indx_sorted,:))
        
    end
end

% if sum(ops.stat.save_est_samp)
%     save([ops.file_dir '\' ops.paradigm_type '_' ops.stat.save_est_samp '.mat'], 'stat_save', 'ops');
% end
   
end