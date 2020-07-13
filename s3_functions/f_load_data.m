function data = f_load_data(ops, app)
    disp('Loading data...');
    

    if ~exist('app', 'var') || isempty('app')
        fw = uifigure;
        new_fig = 1;
    else
        fw = app.UIFigure;
        new_fig = 0;
    end    
    hw = uiprogressdlg(fw,'Title','Loading data...');
    pause(0.05);

    
    for n_cond = 1:numel(ops.conditions_to_analyze)
        cond_name = ops.conditions{ops.conditions_to_analyze(n_cond)};
        num_dsets = numel(ops.file_names.(cond_name));        
        data.(cond_name).num_dsets = num_dsets;
        data.(cond_name).dset_names = cell(num_dsets,1);
        data.(cond_name).c_data = cell(num_dsets,1);
        data.(cond_name).stim_params = cell(num_dsets,1);
        data.(cond_name).num_cells = zeros(num_dsets,1);
        
        % load data
        for n_file = 1:num_dsets
            filename = [cond_name '_' ops.paradigm_type ops.file_names.(cond_name){n_file} '_vectorized_data.mat'];
            data.(cond_name).dset_names{n_file} = filename;
            
            load(filename, 'c_data', 'stim_params');
            data.(cond_name).stim_params{n_file} = stim_params;
            data.(cond_name).num_cells(n_file) = size(c_data.data_vec,1);
            
            % create smooth the spike inference
            % gaussian kernel smoothing
            gaus_win = ops.spike_inf_gaus_win;
            gaus_sig=ops.spike_inf_gaus_sigma; % in frames
            gaus_kernel = exp(-((gaus_win).^2)/(2*gaus_sig^2));
            gaus_kernel = gaus_kernel/sum(gaus_kernel);
            c_data.smooth_spike_inf = conv2(c_data.spike_inf, gaus_kernel, 'same');
            
            data.(cond_name).c_data{n_file} = c_data;
          
            hw.Value = n_file/num_dsets/numel(ops.conditions_to_analyze)  + (n_cond-1)/numel(ops.conditions_to_analyze);
        end
        data.(cond_name).num_cells_total = sum(data.(cond_name).num_cells);
    end

    
    if new_fig
        close(fw);
    else
        close(hw);
    end
end