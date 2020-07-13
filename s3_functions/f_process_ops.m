function ops = f_process_ops(ops)
    ops.conditions = fieldnames(ops.file_names);
    % check if conditions have any datasets
    cond_removed = [];
    for n_cond = 1:numel(ops.conditions_to_analyze)        
        if isempty(ops.file_names.(ops.conditions{ops.conditions_to_analyze(n_cond)}))
            temp_error = [ops.conditions{ops.conditions_to_analyze(n_cond)} '(conditions_to_analyze = ' num2str(ops.conditions_to_analyze(n_cond)) ') is empty and removed drom analysis'];
            disp(['Script error: ' temp_error])
            ops.errors = {ops.errors; temp_error};
            cond_removed = [cond_removed; ops.conditions_to_analyze(n_cond)];
        end
    end
    ops.conditions_to_analyze(cond_removed) = [];
    
    load([ops.conditions{1} '_' ops.paradigm_type ops.file_names.(ops.conditions{1}){1} '_vectorized_data.mat'], 'c_data');
    load_ops = load([ops.conditions{1} '_' ops.paradigm_type ops.file_names.(ops.conditions{1}){1} '_vectorized_data.mat'], 'ops');

    if isfield(load_ops.ops, 'OnACID')
        ops.OnACID = load_ops.ops.OnACID;
    else
        ops.OnACID = 0;
    end
    
    if isfield(load_ops.ops, 'dims')
        ops.dims = load_ops.ops.dims;
    else
        ops.dims = [256, 256];
    end
    
    ops.time_stim_window = c_data.time_stim_window;
    ops.frame_period = c_data.frame_period;
    ops.num_frames = numel(c_data.time_stim_window);
    
    % compute frequencies used
    ops.control_carrier_freq = zeros(1, ops.stim.num_freqs);
    ops.control_carrier_freq(1) = ops.stim.start_freq;
    for ii = 2:ops.stim.num_freqs
        ops.control_carrier_freq(ii) = ops.control_carrier_freq(ii-1) * ops.stim.increase_factor;
    end

    % create freq legend for plots
    ops.context_type_legend = cell(10,3);
    for ii = 1:numel(ops.control_carrier_freq)
        ops.context_type_legend{ii,1} = sprintf('%.1fkHz',ops.control_carrier_freq(ii)/1000);
    end
    ops.context_type_legend{1,2} = 'Cont';
    ops.context_type_legend{1,3} = 'ContFlip';
    ops.context_type_legend{10,2} = 'Dev';
    ops.context_type_legend{10,3} = 'DevFlip';
    for ii = 1:8
        ops.context_type_legend{ii+1,2} = sprintf('Red%d', ii);
        ops.context_type_legend{ii+1,3} = sprintf('Red%dFlip', ii);
    end

    % other parameters
    ops.win.onset_window_frames = logical((ops.time_stim_window>ops.onset_window_time(1)) .* (ops.time_stim_window<ops.onset_window_time(2)));
    ops.win.offset_window_frames = logical((ops.time_stim_window>ops.offset_window_time(1)) .* (ops.time_stim_window<ops.offset_window_time(2)));
    ops.win.baseline_window_frames = logical((ops.time_stim_window>ops.baseline_window_time(1)) .* (ops.time_stim_window<ops.baseline_window_time(2)));
    ops.win.resp_window_frames = logical((ops.time_stim_window>ops.resp_window_time(1)) .* (ops.time_stim_window<ops.resp_window_time(2)));

    % plot 
    ops.win.analysis_window = ops.time_stim_window>ops.baseline_window_time(1);
    ops.analysis_t = ops.time_stim_window(ops.win.analysis_window);


    ops.win.no_base_window = ops.time_stim_window>=0;
    ops.win.no_base_window_size = sum(ops.win.no_base_window);
    
    ops.context_types_all= [1:10, 200+(1:8), 170, 100+(1:8), 270]';        

end