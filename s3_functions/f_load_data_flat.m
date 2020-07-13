function data = f_load_data_flat(ops, app)
    disp('Loading data...');
    num_frames = numel(ops.time_stim_window);
    
    if exist('app', 'var') || isempty('app')
        fw = app.UIFigure;
        new_fig = 0;
    else
        fw = uifigure;
        new_fig = 1;
    end    
    hw = uiprogressdlg(fw,'Title','Loading data...');
    pause(0.05);

    for n_cond = 1:numel(ops.conditions_to_analyze)
        cond_name = ops.conditions{ops.conditions_to_analyze(n_cond)};
        num_dsets = numel(ops.file_names.(cond_name));
        % load data
        % reads the voltage recording file
        temp_data_load = cell(num_dsets,1);
        num_cells_sets = zeros(num_dsets,1);

        num_cells = 0;
        for n_file = 1:num_dsets
            filename = [cond_name '_' ops.paradigm_type ops.file_names.(cond_name){n_file} '_vectorized_data.mat'];
            load(filename, 'c_data', 'stim_params');

            num_cells_sets(n_file) = size(c_data.data_vec,1);

            % here i am taking 42 bins, but ill need to fix this probably
            % by adjusting for different FPS collection speed
            temp_data_load{n_file}.data_vec = c_data.data_vec(:,1:ops.num_frames,:);
            temp_data_load{n_file}.trial_type = c_data.trial_types;
            temp_data_load{n_file}.loco_trials = logical(c_data.loco_trials);
            if isfield(c_data, 'MMN_orientations')
                c_data.MMN_freq = c_data.MMN_orientations;
            elseif isfield(stim_params, 'MMN_freq')
                c_data.MMN_freq = stim_params.MMN_freq;
            end   
            temp_data_load{n_file}.MMN_freq = c_data.MMN_freq;
            temp_data_load{n_file}.A_spac_comp = c_data.A_spac_comp;
            temp_data_load{n_file}.ddt_cell_dataN = c_data.ddt_cell_dataN;
            if ops.OnACID
                temp_data_load{n_file}.data_vec_si = c_data.data_vec_si(:,1:ops.num_frames,:);
                temp_data_load{n_file}.spike_inf = c_data.spike_inf;
            end

            num_cells = num_cells + num_cells_sets(n_file);

            hw.Value = n_file/num_dsets/2/numel(ops.conditions_to_analyze)  + (n_cond-1)/numel(ops.conditions_to_analyze);
        end

        % combine data
        data.(cond_name).vec_data = zeros(num_cells, num_frames, 800);
        data.(cond_name).trial_type = zeros(num_cells, 800);
        data.(cond_name).MMN_freq = zeros(num_cells, 2);
        data.(cond_name).num_cells = num_cells;
        data.(cond_name).num_dataset = zeros(num_cells, 1);
        data.(cond_name).A_spac_comp = zeros(size(temp_data_load{1}.A_spac_comp,1),num_cells);
        data.(cond_name).ddt_cell_dataN = cell(num_cells,1);
        if ops.OnACID
            data.(cond_name).vec_data_si = zeros(num_cells, num_frames, 800);
            data.(cond_name).spike_inf = cell(num_cells,1);
        end
        data.(cond_name).num_cells_sets = num_cells_sets;

        start_cell = 1;

        for n_file = 1:num_dsets
            % get the end cell number
            end_cell = start_cell + size(temp_data_load{n_file}.data_vec,1) - 1;

            data.(cond_name).vec_data(start_cell:end_cell,:,:) = temp_data_load{n_file}.data_vec(:,:,:);
            if ops.OnACID
                data.(cond_name).vec_data_si(start_cell:end_cell,:,:) = temp_data_load{n_file}.data_vec_si(:,:,:);
            end
            data.(cond_name).num_dataset(start_cell:end_cell) = n_file;
            temp_trial_type = temp_data_load{n_file}.trial_type;
            for n_cell = start_cell:end_cell
                data.(cond_name).MMN_freq(n_cell,:) = temp_data_load{n_file}.MMN_freq;
                if ops.remove_loco_trials
                    data.(cond_name).trial_type(n_cell,:) = temp_trial_type.*~temp_data_load{n_file}.loco_trials;
                else
                    data.(cond_name).trial_type(n_cell,:) = temp_trial_type;
                end
            end
            data.(cond_name).A_spac_comp(:,start_cell:end_cell) = temp_data_load{n_file}.A_spac_comp;

            for n_cell = start_cell:end_cell
                data.(cond_name).ddt_cell_dataN{n_cell} = temp_data_load{n_file}.ddt_cell_dataN(n_cell+1-start_cell,:);
                if ops.OnACID
                    data.(cond_name).spike_inf{n_cell} = temp_data_load{n_file}.spike_inf(n_cell+1-start_cell,:);
                end
            end

            start_cell = end_cell + 1;

            hw.Value = (n_file+num_dsets)/num_dsets/2/numel(ops.conditions_to_analyze) + (n_cond-1)/numel(ops.conditions_to_analyze);
        end
    end

    
    if new_fig
        close(fw);
    else
        close(hw);
    end
end