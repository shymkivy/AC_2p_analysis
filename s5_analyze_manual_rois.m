
if 1
    % here need to do stim trig ave
    if sum(strcmpi([AC_data.paradigm], 'behavior'))
        bh_proc = proc_data{strcmpi([AC_data.paradigm], 'behavior')};
        bh_freq_ops = bh_proc.data.stim_params.ops;
        stim_tone_bh = bh_freq_ops.dev_tone_list;
    end
    
    if sum(strcmpi([AC_data.paradigm], 'ammn_stim'))
        ammn_stim_proc = proc_data{strcmpi([AC_data.paradigm], 'ammn_stim')};
        ammn_freq_ops = ammn_stim_proc.data.stim_params.ops;
        stim_run = ammn_freq_ops.stim_trials_volt{1};
        mmn_pat = ammn_freq_ops.MMN_patterns(ammn_freq_ops.paradigm_MMN_pattern(stim_run),:);
        stim_tone_dd = (stim_run - 1)*100 + 70;
        stim_tone_cont = mmn_pat(2);
    end
    
    num_frames_cumsum = cumsum(num_frames(n_pl,:));
    win1 = [5 60];
    
    num_paradigms = numel(AC_data.paradigm);
    for n_pr = 1:num_paradigms
        if sum(strcmpi(AC_data.paradigm{n_pr}, {'behavior', 'ammn'}))
            stim_chan = 'stim type';
        elseif sum(strcmpi(AC_data.paradigm{n_pr}, {'ammn_stim', 'spont_stim'}))
            stim_chan = 'pockel';
        elseif strcmpi(AC_data.paradigm{n_pr}, 'spont')
            stim_chan = '';
        else
            fprintf('paradigm %s not listed\n', AC_data.paradigm{n_pr});
            stim_chan = '';
        end
        
        if numel(stim_chan)
            idx_chan = strcmpi(proc_data{n_pr}.ops.chan_labels, stim_chan);
            stim_times_frame = proc_data{n_pr}.data.stim_times_frame{idx_chan,n_pl};
            
            stim_times_frame(stim_times_frame < win1(1)) = [];
            stim_times_frame((stim_times_frame + win1(2)) > num_frames(n_pr)) = [];
            
            stim_trace2 = trace1n((num_frames_cumsum(n_pr)-num_frames(n_pl,n_pr)+1):num_frames_cumsum(n_pr));
            sorted_trials = squeeze(f_get_stim_trig_resp(stim_trace2', stim_times_frame, win1));
            stim_trace2d = trace1nd((num_frames_cumsum(n_pr)-num_frames(n_pl,n_pr)+1):num_frames_cumsum(n_pr));
            sorted_trialsd = squeeze(f_get_stim_trig_resp(stim_trace2d, stim_times_frame, win1));

            figure; 
            subplot(3,1,1);
            imagesc(sorted_trials'); 
            title('trials all');
            subplot(3,1,2); hold on;
            plot(sorted_trials, 'color', [.6 .6 .6])
            plot(mean(sorted_trials,2), 'Linewidth', 2, 'color', 'k');
            axis tight;
            title('trig ave');
            subplot(3,1,3); hold on;
            plot(sorted_trialsd, 'color', [.6 .6 .6])
            plot(mean(sorted_trialsd,2), 'Linewidth', 2, 'color', 'k');
            axis tight;
            title('trig ave deconvolved');
            sgtitle(stim_chan);
        end
        
        if strcmpi(AC_data.paradigm{n_pr}, 'ammn')
            1
            
        end
        
    end
    
    
    
    idx_stim = strcmpi([AC_data.paradigm], 'ammn_stim');
    if sum(idx_stim)
        
    end
    
    idx_stim = strcmpi([AC_data.paradigm], 'behavior');
    if sum(idx_stim)
        idx_chan = strcmpi(proc_data{idx_stim}.ops.chan_labels, 'stim type');
        stim_times_frame = proc_data{idx_stim}.data.stim_times_frame{idx_chan,n_pl};
        
        stim_trace2 = trace1n((num_frames_cumsum(idx_stim)-num_frames(n_pl,idx_stim)+1):num_frames_cumsum(idx_stim));
        sorted_trials = squeeze(f_get_stim_trig_resp(stim_trace2', stim_times_frame, win1));
        
        figure; 
        subplot(2,1,1);
        imagesc(sorted_trials');
        title('trials all');
        subplot(2,1,2); hold on;
        plot(sorted_trials, 'color', [.6 .6 .6])
        plot(mean(sorted_trials,2), 'Linewidth', 2, 'color', 'k');
        axis tight;
        title('trig ave');
        sgtitle('behavior');
    end
    
    idx_stim = strcmpi([AC_data.paradigm], 'ammn');
    if sum(idx_stim)
        idx_stim = strcmpi(proc_data{idx_stim}.ops.chan_labels, 'stim type');
        stim_times_frame = proc_data{idx_stim}.data.stim_times_frame{idx_chan,n_pl};

    end
    % analyze tuning before + a
    
end