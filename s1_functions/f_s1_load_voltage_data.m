function data = f_s1_load_voltage_data(data, ops)
frame_data_XML = data.frame_data_XML;
volt_data_all_segmented = cell(ops.num_volt_in_files,1);
volt_times_all_segmented = cell(ops.num_volt_in_files,1);

% load and process voltage files
for n_file = 1:ops.num_volt_in_files
    dattmp = dlmread([ops.file_dir '\' ops.files_volt_in{n_file} '.csv'],',',1,0);
    
    volt_time = dattmp(:,1);
    dattmp2 = dattmp(:,2:end);
    
    [T, num_volt_chan] = size(dattmp2);
        
    % syncs the start of voltage recording with start of frame
    %~strcmp(data.frame_data_XML.trigger_mode{n_file}, 'Start with next scan (PFI0)')
    if ~sum(strfind(data.frame_data_XML.trigger_mode{n_file}, 'next scan'), strfind(data.frame_data_XML.trigger_mode{n_file}, 'external trigger'))
        dattmp2 = dattmp2(round((frame_data_XML.frame_start(n_file)-frame_data_XML.volt_start(n_file))+1):end,:);
    end

    % reorder chans
    all_chan = 1:num_volt_chan;
    idx1 = logical(sum(all_chan == ops.volt_chan_order'));
    idx2 = logical(sum(all_chan == ops.volt_chan_order',2)); 
    chan_list = [ops.volt_chan_order(idx2) all_chan(~idx1)];
    
    % process voltage stim trace
    dat_proc = dattmp2(:,chan_list);
   
    med_filt_list = [1 2 5 6];
    remove_base_list = [1 2 4 5 6];
    cap_lim = [2];
    
    for n_ch = 1:num_volt_chan
        temp_trace = dat_proc(:, n_ch);
        if sum(n_ch == med_filt_list)
            temp_trace = medfilt1(temp_trace, 49);
        end
        if sum(n_ch == remove_base_list)
            temp_trace = temp_trace - mode(temp_trace);
        end
        if sum(n_ch == cap_lim)
            temp_trace(temp_trace>1) = 1;
        end
        dat_proc(:, n_ch) = temp_trace;
    end
    
    % make sure the median filter window is odd
%     dat_proc(:,1) = dat_proc(:,1) - median(dat_proc(:,1));
%     dat_proc(:,1) = medfilt1(dat_proc(:,1), 49);
%     
%     dat_proc(:,2) = dat_proc(:,2) - median(dat_proc(:,2));
%     dat_proc(:,2) = medfilt1(dat_proc(:,2), 49);
%     dat_proc(:,2) = if_rescale(dat_proc(:,2));
%     if num_volt_chan > 3
%         dat_proc(:,4) = medfilt1(dat_proc2(:,4), 49);
%     end  
    
    % fill in or cut voltage to make make length the same with XML frames length
    volt_t = mean(diff(volt_time));
    if round(frame_data_XML.frame_times_raw{n_file}(end)) > T
        add_times = round(frame_data_XML.frame_times_raw{n_file}(end) - T);
        dat_proc_pad = [dat_proc; zeros(add_times, num_volt_chan)];
        volt_time_pad = [volt_time; (volt_time(end)+(1:add_times)*volt_t)'];
    else
        dat_proc_pad = dat_proc(1:round(frame_data_XML.frame_times_raw{n_file}(end)),:);
        volt_time_pad = volt_time(1:round(frame_data_XML.frame_times_raw{n_file}(end)));
    end
    
    
    volt_data_all_segmented{n_file} = dat_proc_pad;
    volt_times_all_segmented{n_file} = volt_time_pad;
end
data.volt_data_all_segmented = volt_data_all_segmented;
data.volt_data_all = cat(1,volt_data_all_segmented{:});
data.num_samp_volt = size(data.volt_data_all,1);
data.volt_time_all = cat(1,volt_times_all_segmented{:});



if ops.plot_details
    figure;
    hold on;
    sp1 = cell(num_volt_chan,1);
    for n_plt = 1:num_volt_chan
        sp1{n_plt} = subplot(num_volt_chan,1,n_plt);
        plot(dat_proc(:,n_plt));
        if n_plt <= numel(ops.volt_chan_labels)
            title(ops.volt_chan_labels{n_plt})
        end
        axis tight;
    end
    sgtitle(['Using channel ' num2str(ops.align_to_channel) ' for alignment']);
    drawnow();
    linkaxes(cat(1,sp1{:}), 'x')
end

end

function trace_out = if_rescale(trace)

base = min(trace);
base_sub = trace - base;
peak = max(base_sub);
trace_out = base_sub/peak;

end