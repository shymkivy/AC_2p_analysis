function volt_trace_out = align_volt_trace(volt_trace_in, synch_points)
    % synch_points needs Ca points in column 1 and voltage points in column 2

    % align traces
    % compute shift
    synch_points = round(synch_points);
    volt_delay = synch_points(2,1) - synch_points(1,1);
    if volt_delay > 0 % shift right
        shifted_voltage_trace = [volt_trace_in(volt_delay+1:end); zeros(volt_delay, 1)];
    elseif volt_delay < 0 %  shift left
        shifted_voltage_trace = [zeros(abs(volt_delay), 1); volt_trace_in(1:end-(abs(volt_delay)))];
    else
        shifted_voltage_trace = volt_trace_in;
    end


    % compute change in size
    volt_trace_fix = shifted_voltage_trace;
    for ii = 2:(size(synch_points,2))
        size_change = (synch_points(2,ii) - synch_points(2,ii-1)) - (synch_points(1,ii) - synch_points(1,ii-1));
        add_remove_frames2 = round(linspace(synch_points(1, ii-1), synch_points(1, ii), abs(size_change)+2));
        add_remove_frames = add_remove_frames2(2:end-1);
        if size_change > 0 % to contract voltage trace
            temp_voltage = volt_trace_fix;
            temp_voltage2 = temp_voltage;
            temp_voltage2(add_remove_frames) = [];
            temp_voltage = [temp_voltage2; zeros(length(add_remove_frames), 1)];
            volt_trace_fix = temp_voltage;
        elseif size_change < 0 % to expand voltage trace
            temp_voltage = zeros([length(volt_trace_fix)+ abs(size_change), 1]);
            insertion_locations = temp_voltage;
            insertion_locations(add_remove_frames) = 1;
            copied_locations = temp_voltage;
            copied_locations(add_remove_frames+1) = 1;
            temp_voltage(logical(1-insertion_locations)) = volt_trace_fix;
            temp_voltage(logical(insertion_locations)) = temp_voltage(logical(copied_locations));
            volt_trace_fix = temp_voltage(1:length(volt_trace_fix));
        end
    end
    volt_trace_out = volt_trace_fix;
end