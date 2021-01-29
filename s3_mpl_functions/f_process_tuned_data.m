function data = f_process_tuned_data(data, ops)

num_dset = numel(data.area);

for n_dset = 1:num_dset   
    %% full
    tuning_temp = data.tuning_all{n_dset}.peak_tuning_full_resp;
    fr_peak_tuned_trials_full = tuning_temp.fr_peak_mag_ave_z>=ops.stat.z_scores_thresh;
    fr_peak_tuned_trials_full = fr_peak_tuned_trials_full.*(tuning_temp.fr_peak_reliability>=ops.stat.reliability_thresh);

    data.peak_tuned_trials_full{n_dset,1} = fr_peak_tuned_trials_full;
    data.peak_tuned_trials_full_ctx{n_dset,1} = fr_peak_tuned_trials_full(:,data.ctx_mmn{n_dset});
    data.peak_tuned_trials_full_reliab{n_dset,1} = data.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_reliability;

    %% onset
    tuning_temp = data.tuning_all{n_dset}.peak_tuning_onset;
    fr_peak_tuned_trials_onset = tuning_temp.fr_peak_mag_ave_z>=ops.stat.z_scores_thresh;
    fr_peak_tuned_trials_onset = fr_peak_tuned_trials_onset.*(tuning_temp.fr_peak_reliability>=ops.stat.reliability_thresh);

    data.peak_tuned_trials_onset{n_dset,1} = fr_peak_tuned_trials_onset;
    data.peak_tuned_trials_onset_ctx{n_dset,1} = fr_peak_tuned_trials_onset(:,data.ctx_mmn{n_dset});
    data.peak_tuned_trials_onset_reliab{n_dset,1} = data.tuning_all{n_dset}.peak_tuning_onset.fr_peak_reliability;

    %% offset
    tuning_temp = data.tuning_all{n_dset}.peak_tuning_offset;
    fr_peak_tuned_trials_offset = tuning_temp.fr_peak_mag_ave_z>=ops.stat.z_scores_thresh;
    fr_peak_tuned_trials_offset = fr_peak_tuned_trials_offset.*(tuning_temp.fr_peak_reliability>=ops.stat.reliability_thresh);

    data.peak_tuned_trials_offset{n_dset,1} = fr_peak_tuned_trials_offset;
    data.peak_tuned_trials_offset_ctx{n_dset,1} = fr_peak_tuned_trials_offset(:,data.ctx_mmn{n_dset});
    data.peak_tuned_trials_offset_reliab{n_dset,1} = data.tuning_all{n_dset}.peak_tuning_offset.fr_peak_reliability;

    %% combined
    fr_pk_tuned_trials_combined = logical(fr_peak_tuned_trials_full+fr_peak_tuned_trials_onset+fr_peak_tuned_trials_offset);

    data.peak_tuned_trials_combined{n_dset,1} = fr_pk_tuned_trials_combined;
    data.peak_tuned_trials_combined_ctx{n_dset,1} = fr_pk_tuned_trials_combined(:,data.ctx_mmn{n_dset});


end

end