function [data, ops] = f_s1_load_voltage_data(data, ops)
frame_data_XML = data.frame_data_XML;
volt_data_all_segmented = cell(ops.num_volt_in_files,1);
volt_times_all_segmented = cell(ops.num_volt_in_files,1);

% load and process voltage files
for n_file = 1:ops.num_volt_in_files
    dattmp = dlmread([ops.file_dir '\' ops.files_volt_in{n_file} '.csv'],',',1,0);
    
    volt_time = dattmp(:,1);
    dattmp2 = dattmp(:,2:end);
    
    [T, num_volt_rec] = size(dattmp2);

    % syncs the start of voltage recording with start of frame
    %~strcmp(data.frame_data_XML.trigger_mode{n_file}, 'Start with next scan (PFI0)')
    if ~sum(strfind(data.frame_data_XML.trigger_mode{n_file}, 'next scan'), strfind(data.frame_data_XML.trigger_mode{n_file}, 'external trigger'))
        dattmp2 = dattmp2(round((frame_data_XML.frame_start(n_file)-frame_data_XML.volt_start(n_file))+1):end,:);
    end
    
    idx1 = logical(ops.volt_chan_order);
    ops.chan_list = ops.volt_chan_order(idx1);
    ops.chan_type = repmat({'base_volt'}, sum(idx1), 1);
    ops.chan_idx = ones(sum(idx1),1);
    ops.chan_labels = ops.volt_chan_labels(idx1);
    if isfield(ops, 'volt_chan_order_bh')
        idx1 = logical(ops.volt_chan_order_bh);
        ops.chan_list = [ops.chan_list, ops.volt_chan_order_bh(idx1)];
        ops.chan_type = [ops.chan_type; repmat({'behavior'},sum(idx1), 1)];
        ops.chan_idx = [ops.chan_idx; ones(sum(idx1),1)*2];
        ops.chan_labels = [ops.chan_labels, ops.volt_chan_labels_bh(idx1)];
    end
    if isfield(ops, 'volt_chan_order_stim')
        idx1 = logical(ops.volt_chan_order_stim);
        ops.chan_list = [ops.chan_list, ops.volt_chan_order_stim(idx1)];
        ops.chan_type = [ops.chan_type; repmat({'opto_stim'},sum(idx1), 1)];
        ops.chan_idx = [ops.chan_idx; ones(sum(idx1),1)*3];
        ops.chan_labels = [ops.chan_labels, ops.volt_chan_labels_stim(idx1)];
    end
    
%     % reorder chans
%     all_chan = 1:num_volt_rec;
%     idx1 = logical(sum(all_chan == ops.volt_chan_order'));
%     idx2 = logical(sum(all_chan == ops.volt_chan_order',2)); 
%     chan_list = [ops.volt_chan_order(idx2) all_chan(~idx1)];
%     
    % process voltage stim trace
    dat_proc = dattmp2(:,ops.chan_list);
    
    % 1 stim
    % 2 led
    % 3 loco
    % 4 auditory TDT
    % 5 lick
    % 6 reward
    % 7 LED bh
    % 8 pockel
    % 9 slm_stim_type
    
    med_filt_list = [1 2 5 6 7 8 9]; 
    remove_base_list = [1 2 4 5 6 7 8 9];
    cap_lim = 2;
    
    for n_ch = 1:numel(ops.chan_list)
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
    
    unused_chan = 1:num_volt_rec;
    unused_chan(ops.chan_list) = [];
    
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
        dat_proc_pad = [dat_proc; zeros(add_times, numel(ops.chan_list))];
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
    if numel(unused_chan)
        figure;
        hold on;
        for n_plt = 1:numel(unused_chan)
            subplot(numel(unused_chan),1,n_plt);
            plot(dattmp2(:,unused_chan(n_plt)));
            title(sprintf('chan %d', unused_chan(n_plt)));
        end
        sgtitle(sprintf('Unused channels; %s', ops.file_core), 'interpreter', 'none');
    end
    
    figure;
    hold on;
    sp1 = cell(numel(ops.chan_list),1);
    for n_plt = 1:numel(ops.chan_list)
        sp1{n_plt} = subplot(numel(ops.chan_list),1,n_plt);
        plot(dat_proc(:,n_plt));
        title(sprintf('%s; %s; chan %d', ops.chan_labels{n_plt}, ops.chan_type{n_plt}, ops.chan_list(n_plt)), 'interpreter', 'none')
        axis tight;
    end
    sgtitle(sprintf('Using channel %d for alignment; %s', ops.align_to_channel, ops.file_core), 'interpreter', 'none');
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