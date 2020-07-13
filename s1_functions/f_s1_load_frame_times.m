function [data, ops] = f_s1_load_frame_times(data, ops)
% each imaged file has its own xml with frame times. If multiplane imaging
% was used, those need to be split up into the each plane. if multiple xml
% files were input, the ca input is concatenated across them
%
%   variables(num_volt_files, num_planes)
%

%% XML first
if ~isfield(data, 'frame_data_XML')    
    % initialize parameters, yesss
    frame_data_XML.frame_period_raw = zeros(ops.num_volt_in_files,1);
    frame_data_XML.volt_start = zeros(ops.num_volt_in_files,1);
    frame_data_XML.frame_start = zeros(ops.num_volt_in_files,1);
    frame_data_XML.num_frames_raw = zeros(ops.num_volt_in_files,1);
    frame_data_XML.frame_times_raw = {ops.num_volt_in_files,1};
    frame_data_XML.trigger_mode = {ops.num_volt_in_files,1};
    for n_file = 1:ops.num_volt_in_files
        fprintf('Loading file %d of %d...\n', n_file, ops.num_volt_in_files);
        % use frame times from XML file
        
        % extract frame data from XML file with function, obviously
        data_out = extract_frame_data_from_XML2([ops.file_dir '\' ops.files_volt_in{n_file}]);

        % move data
        frame_data_XML.frame_period_raw(n_file) = data_out.frame_period;
        frame_data_XML.volt_start(n_file) = data_out.volt_start;
        frame_data_XML.frame_start(n_file) = data_out.frame_start;
        frame_data_XML.num_frames_raw(n_file) = data_out.total_frames;
        % the signal is captured at end of frame, so shift everything by period (11/25/19)
        frame_data_XML.frame_times_raw{n_file} = data_out.frame_times+data_out.frame_period;
        frame_data_XML.trigger_mode{n_file} = data_out.trigger_mode;
        
        if ~isempty(data_out.errors)
            ops.errors = [ops.errors; data_out.errors];
        end
    end
    data.frame_data_XML = frame_data_XML;
    if exist(ops.file_save_path_full_processing_params, 'file')
        save(ops.file_save_path_full_processing_params, 'frame_data_XML', '-append')
    else
        save(ops.file_save_path_full_processing_params, 'frame_data_XML')
    end
    fprintf('Saved frame_data\n');
else
    frame_data_XML = data.frame_data_XML;
end


%% process the frame data
frame_data.num_frames_mpl = zeros(ops.num_volt_in_files,ops.num_planes);
frame_data.frame_times_mpl = cell(ops.num_volt_in_files,ops.num_planes);
frame_data.volume_period = zeros(ops.num_volt_in_files,1);
for n_file = 1:ops.num_volt_in_files
    frame_data.volume_period(n_file) = frame_data_XML.frame_period_raw(n_file)*ops.num_planes;
    [frame_data.frame_times_mpl(n_file,:), frame_data.num_frames_mpl(n_file,:)] = f_s1_multiplane_split(frame_data_XML.frame_times_raw{n_file}, ops.num_planes);
end
frame_data.volume_period_ave = mean(frame_data.volume_period(:));
frame_data.frame_period_ave = mean(frame_data_XML.frame_period_raw(:));
frame_data.num_frames_all = sum(frame_data.num_frames_mpl(:));

data.frame_data = frame_data;
end