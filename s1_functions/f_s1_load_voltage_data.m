function data = f_s1_load_voltage_data(data, ops)
frame_data_XML = data.frame_data_XML;
volt_data_all_segmented = cell(ops.num_volt_in_files,1);

% load and process voltage files
for n_file = 1:ops.num_volt_in_files
    dattmp = dlmread([ops.file_dir '\' ops.files_volt_in{n_file} '.csv'],',',1,1);
    num_volt_traces = size(dattmp,2);
    
    
    % syncs the start of voltage recording with start of frame
    %~strcmp(data.frame_data_XML.trigger_mode{n_file}, 'Start with next scan (PFI0)')
    if ~sum(strfind(data.frame_data_XML.trigger_mode{n_file}, 'next scan'), strfind(data.frame_data_XML.trigger_mode{n_file}, 'external trigger'))
        dattmp = dattmp(round((frame_data_XML.frame_start(n_file)-frame_data_XML.volt_start(n_file))+1):end,:);
    end
    
    % fill in or cut voltage to make make length the same with XML frames length
    num_pts = size(dattmp,1);
    if round(frame_data_XML.frame_times_raw{n_file}(end)) > num_pts
        dattmp_pad = [dattmp; zeros(round(frame_data_XML.frame_times_raw{n_file}(end) - num_pts), size(dattmp,2))];
    else
        dattmp_pad = dattmp(1:round(frame_data_XML.frame_times_raw{n_file}(end)),:);
    end

    % process voltage stim trace
    dat_proc = zeros(size(dattmp_pad,1),num_volt_traces);
    % make sure the median filter window is odd
    dat_proc(:,1) = medfilt1(dattmp_pad(:,ops.parameters.stimchan), 49);
    temp_trace = medfilt1(dattmp_pad(:,ops.parameters.ledchan), 49);
    dat_proc(:,2) = if_normalize(temp_trace);
    dat_proc(:,3) = dattmp_pad(:,ops.parameters.movchan);
    if num_volt_traces > 3
        dat_proc(:,4) = medfilt1(dattmp_pad(:,ops.parameters.TDT_volt_chan), 49);
    end  
    
    volt_data_all_segmented{n_file} = dat_proc;
end
data.volt_data_all_segmented = volt_data_all_segmented;
data.volt_data_all = cat(1,volt_data_all_segmented{:});
data.num_samp_volt = size(data.volt_data_all,1);



if ops.plot_details
    figure;
    hold on;
    for n_plt = 1:num_volt_traces
        subplot(num_volt_traces,1,n_plt);
        plot(dattmp(:,n_plt));
        title(ops.volt_chan_labels{n_plt})
    end
    sgtitle(['Using channel ' num2str(ops.align_to_channel) ' for alignment']);
    drawnow();
end

end

function norm_trace = if_normalize(trace)
base = min(trace);
base_sub = trace - base;
peak = max(base_sub);
norm_trace = base_sub/peak;
end