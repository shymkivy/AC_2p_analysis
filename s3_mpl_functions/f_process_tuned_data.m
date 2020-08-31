function data = f_process_tuned_data(data, ops)


for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
%     if ~load_file
%         %stat_save.(cond_name).onset_thresh = cell(data.(cond_name).num_dsets,1);
%         %stat_save.(cond_name).offset = cell(data.(cond_name).num_dsets,1);
%     end
    for n_dset = 1:cdata.num_dsets       
        %% full
        tuning_temp = cdata.tuning_all{n_dset}.peak_tuning_full_resp;
        fr_peak_tuned_trials_full = tuning_temp.fr_peak_mag_ave_z>=ops.stat.z_scores_thresh;
        fr_peak_tuned_trials_full = fr_peak_tuned_trials_full.*(tuning_temp.fr_peak_reliability>=ops.stat.reliability_thresh);
        
        cdata.peak_tuned_trials_full{n_dset,1} = fr_peak_tuned_trials_full;
        cdata.peak_tuned_trials_full_ctx{n_dset,1} = fr_peak_tuned_trials_full(:,cdata.ctx_mmn{n_dset});
        cdata.peak_tuned_trials_full_reliab{n_dset,1} = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_reliability;

        %% onset
        tuning_temp = cdata.tuning_all{n_dset}.peak_tuning_onset;
        fr_peak_tuned_trials_onset = tuning_temp.fr_peak_mag_ave_z>=ops.stat.z_scores_thresh;
        fr_peak_tuned_trials_onset = fr_peak_tuned_trials_onset.*(tuning_temp.fr_peak_reliability>=ops.stat.reliability_thresh);
        
        cdata.peak_tuned_trials_onset{n_dset,1} = fr_peak_tuned_trials_onset;
        cdata.peak_tuned_trials_onset_ctx{n_dset,1} = fr_peak_tuned_trials_onset(:,cdata.ctx_mmn{n_dset});
        cdata.peak_tuned_trials_onset_reliab{n_dset,1} = cdata.tuning_all{n_dset}.peak_tuning_onset.fr_peak_reliability;
        
        %% offset
        tuning_temp = cdata.tuning_all{n_dset}.peak_tuning_offset;
        fr_peak_tuned_trials_offset = tuning_temp.fr_peak_mag_ave_z>=ops.stat.z_scores_thresh;
        fr_peak_tuned_trials_offset = fr_peak_tuned_trials_offset.*(tuning_temp.fr_peak_reliability>=ops.stat.reliability_thresh);
        
        cdata.peak_tuned_trials_offset{n_dset,1} = fr_peak_tuned_trials_offset;
        cdata.peak_tuned_trials_offset_ctx{n_dset,1} = fr_peak_tuned_trials_offset(:,cdata.ctx_mmn{n_dset});
        cdata.peak_tuned_trials_offset_reliab{n_dset,1} = cdata.tuning_all{n_dset}.peak_tuning_offset.fr_peak_reliability;

        %% combined
        fr_pk_tuned_trials_combined = logical(fr_peak_tuned_trials_full+fr_peak_tuned_trials_onset+fr_peak_tuned_trials_offset);
        
        cdata.peak_tuned_trials_combined{n_dset,1} = fr_pk_tuned_trials_combined;
        cdata.peak_tuned_trials_combined_ctx{n_dset,1} = fr_pk_tuned_trials_combined(:,cdata.ctx_mmn{n_dset});
        

    end
    data.(cond_name) = cdata;
end







end