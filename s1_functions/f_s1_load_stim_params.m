function [data, ops] = f_s1_load_stim_params(data, ops)
    
    if exist([ops.file_dir '\' ops.file_core '_stim_data.mat'], 'file')
    %     load([file '_stim_data.mat'], 'cont_trials_seq', 'MMN_trials_seq', 'MMN_orientations');
        stim_params = load([ops.file_dir '\' ops.file_core '_stim_data.mat']);
        
        if isfield(stim_params, 'ops')
            if isfield(stim_params.ops, 'stim_time')
                stim_params.stim_duration = stim_params.ops.stim_time;
            end
            if isfield(stim_params.ops, 'isi_time')
                stim_params.isi = stim_params.ops.isi_time;
            end
            if isfield(stim_params.ops, 'num_freqs')
                stim_params.num_freqs = stim_params.ops.num_freqs;
            end
            if isfield(stim_params.ops, 'start_freq')
                stim_params.start_freq = stim_params.ops.start_freq;
            end
            if isfield(stim_params.ops, 'end_freq')
                stim_params.end_freq = stim_params.ops.end_freq;
            end
            if isfield(stim_params.ops, 'increase_factor')
                stim_params.increase_factor = stim_params.ops.increase_factor;
            end
            if isfield(stim_params.ops, 'paradigm_sequence')
                data.MMN_orientations = stim_params.ops.MMN_patterns(stim_params.ops.paradigm_MMN_pattern(strcmpi(stim_params.ops.paradigm_sequence,'mmn')),:);
            end
        end
        
        if isfield(stim_params, 'duration')
            stim_params.stim_duration = stim_params.duration;
            stim_params = rmfield(stim_params, 'duration');
        end
        
        if isfield(stim_params, 'freq_grating_stim')
            stim_params = rmfield(stim_params, 'freq_grating_stim');
        end
        
        if isfield(stim_params, 'MMN_orientations')
            data.MMN_orientations = stim_params.MMN_orientations;
            stim_params = rmfield(stim_params, 'MMN_orientations');
        end
        
        if isfield(stim_params, 'MMN_trials_seq')
            data.MMN_trials_seq = stim_params.MMN_trials_seq;
            stim_params = rmfield(stim_params, 'MMN_trials_seq');
        end
        
        if isfield(stim_params, 'cont_trials_seq')
            data.cont_trials_seq = stim_params.cont_trials_seq;
            stim_params = rmfield(stim_params, 'cont_trials_seq');
        end
        
        if isfield(stim_params, 'freq_grating_stim_dsp')
            data.freq_grating_stim_dsp = stim_params.freq_grating_stim_dsp;
            stim_params = rmfield(stim_params, 'freq_grating_stim_dsp');
        end
        
        if ~isfield(stim_params, 'stim_duration')
            stim_params.stim_duration = 0.5;
            ops.errors = [ops.errors; 'manual written stim_params.stim_duration = ' num2str(stim_params.stim_duration)];
        end
        
        if ~isfield(stim_params, 'isi')
            stim_params.isi = 0.5;
            ops.errors = [ops.errors; 'manual written stim_params.isi = ' num2str(stim_params.isi)];
        end
        
        if isfield(stim_params, 'grating_angles')
            stim_params.num_freqs = numel(stim_params.grating_angles);
        end
        
        if ~isfield(stim_params, 'num_freqs')
            if isfield(stim_params.ops, 'angs_rad')
                stim_params.num_freqs = numel(stim_params.ops.angs_rad);
            else
                error('Extract num freqs from stim params')
            end
        end
    else
    %     stim_params.MMN_orientations = [4, 7];
        stim_params.stim_duration = 0.5;
        stim_params.isi = 0.5;
        stim_params.start_freq = 2000; %kHz
        stim_params.increase_factor = 1.5;
        stim_params.num_freqs = 10;
        ops.errors = [ops.errors; 'manual written stim_params'];
    end
    data.stim_params = stim_params;
end
