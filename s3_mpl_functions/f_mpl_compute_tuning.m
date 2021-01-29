function data = f_mpl_compute_tuning(data, ops)

disp('Computing tunning...');

if ops.waitbar
    wb = f_waitbar_initialize([], 'Computing tuning...');
end

num_dset = numel(data.area);

for n_dset = 1:num_dset

    trial_data_sort = cell(1,data.num_planes(n_dset));
    trial_data_sort_sm = cell(1,data.num_planes(n_dset));
    trial_num_baseline_resp_frames = data.trial_window{n_dset}.trial_num_baseline_resp_frames;
    for n_pl = 1:data.num_planes(n_dset)
        stim_frame_index = data.stim_frame_index{n_dset,n_pl};
        firing_rate = data.firing_rate{n_dset,n_pl};
        firing_rate_sm = data.firing_rate_smooth{n_dset,n_pl};
        trial_data_sort{n_pl} = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_num_baseline_resp_frames);
        trial_data_sort_sm{n_pl} = f_get_stim_trig_resp(firing_rate_sm, stim_frame_index, trial_num_baseline_resp_frames); 
    end
    trial_data_sort = cat(1,trial_data_sort{:});
    trial_data_sort_sm = cat(1,trial_data_sort_sm{:});
    trial_types = data.trial_types{n_dset};

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
        data.trial_types{n_dset} = trial_types;
    end
    %%
    MMN_freq = data.MMN_freq{n_dset};
    ctx_mmn = [MMN_freq(2), 19, 20, ...
               MMN_freq(1), 29, 30];  

    

    %%
    [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, MMN_freq, ops);
    [trial_data_sort_sm_wctx, ~] =  f_s3_add_ctx_trials(trial_data_sort_sm, trial_types, MMN_freq, ops);
    data.trial_data_sort_wctx{n_dset} = trial_data_sort_wctx;
    data.trial_data_sort_sm_wctx{n_dset} = trial_data_sort_sm_wctx;
    data.trial_types_wctx{n_dset} = trial_types_wctx;

    
    dset_params = data(n_dset,:);
    dset_params.ctx_mmn{1} = ctx_mmn;
    %% get trial averages compute tunning     
    tuning_all = f_get_tuning_dset(trial_data_sort_sm_wctx, ops.context_types_all, dset_params, ops);
    data.tuning_all{n_dset,1} = tuning_all;

    %%
    z_thresh_all = max(nanmean(tuning_all.trace_tuning.stat_trace.z_factors,3),[],2);

    trial_ave = f_mpl_trial_average(trial_data_sort_sm_wctx,trial_types_wctx, ops.context_types_all, 'none');
    trial_ave_z = trial_ave./z_thresh_all;
    %trial_ave = (trial_ave-means_all)./z_thresh_all;

    data.trial_ave{n_dset,1} = trial_ave;
    data.trial_ave_z{n_dset,1} = trial_ave_z;
    data.trial_ave_mmn{n_dset,1} = trial_ave(:,:,ctx_mmn);
    data.trial_ave_mmn_z{n_dset,1} = trial_ave_z(:,:,ctx_mmn);
    %cdata.resp_cells_mmn_onset{n_dset} = data.(cond_name).resp_cells_all_onset{n_dset}(:,ctx_mmn);
    %cdata.resp_cells_mmn_offset{n_dset} = data.(cond_name).resp_cells_all_offset{n_dset}(:,ctx_mmn);
    data.ctx_mmn{n_dset,1} = ctx_mmn;

    if ops.waitbar
        f_waitbar_update(wb, n_dset/num_dset, sprintf('Computing tuning dset %d/%d', n_dset, num_dset));
    end

end

if ops.waitbar
    f_waitbar_close(wb);
end

end