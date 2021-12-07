function [params] = f_s0_compute_align_cuts(params)
    if ~isfield(params, 'auto_align_pulse_crop')
        params.auto_align_pulse_crop = 1; % default
    end
    ave_trace = params.ave_trace;

    min_trace = min(ave_trace);
    max_trace = max(ave_trace);
    norm_ave_trace = (ave_trace - min_trace)/(max_trace-min_trace);
    
    if params.auto_align_pulse_crop
        thresh = 0.5;
        pulse_buff = 30; % frames
        
        pulse_trace = (norm_ave_trace);
        pulse_trace(pulse_trace<thresh) = 0;
        pulse_trace(pulse_trace>thresh) = 1;
       
        pulse_on = find(diff(pulse_trace)>0)+1;
        pulse_off = find(diff(pulse_trace)<0)+1;
        
        % quality check if light turns off in beginning 
        if and(pulse_trace(1) == 1,numel(pulse_off)>numel(pulse_on))
            pulse_off(1) = [];
        end
        
        num_frag = round(numel(pulse_on)-1);
        
        vid_cuts = zeros(num_frag,2);
        
        n_pulse = 1;
        vid_cuts_trace = false(size(norm_ave_trace));
        for n_frag = 1:num_frag
            vid_cuts(n_frag,1) = pulse_off(n_pulse) + pulse_buff;
            n_pulse = n_pulse+1;
            vid_cuts(n_frag,2) = pulse_on(n_pulse) - pulse_buff;
            vid_cuts_trace(vid_cuts(n_frag,1):vid_cuts(n_frag,2)) = 1;
        end

        
        figure; plot(norm_ave_trace);
        hold on; plot(vid_cuts_trace);
        axis tight;
        title(sprintf('Automatic %d frag selected %s', num_frag, params.title_tag), 'interpreter', 'none');
    else
        f1 = figure;
        plot(norm_ave_trace);
        axis tight;
        title('how many fragments?');
        num_frag = input('how many fragments? (int):');
        vid_cuts = zeros(num_frag,2);
        vid_cuts_trace = false(size(norm_ave_trace));
        for n_frag = 1:num_frag
            title(sprintf('Select fragment %d/%d (2 clicks)', n_frag,num_frag));
            [temp_cuts, ~] = ginput(2);
            vid_cuts(n_frag,:) = round(temp_cuts);
            if vid_cuts(n_frag,1) < 1
                vid_cuts(n_frag,1) = 1;
            end
            if vid_cuts(n_frag,2) > numel(ave_trace)
                vid_cuts(n_frag,2) = numel(ave_trace);
            end
            vid_cuts_trace(vid_cuts(n_frag,1):vid_cuts(n_frag,2)) = 1;
            plot(norm_ave_trace);
            hold on;
            plot(vid_cuts_trace);
            hold off;
            axis tight;
        end
        close(f1)
        
        figure;plot(norm_ave_trace);
        hold on; plot(vid_cuts_trace);
        axis tight;
        title(sprintf('Manual %d frag selected %s', num_frag, params.mpl_tag), 'interpreter', 'none');
    end
    
    params.vid_cuts = vid_cuts;
    params.vid_cuts_trace = vid_cuts_trace;
end
