function ops = f_s1_process_ops(ops)

if strcmp(ops.ca_processing ,'onacid')
    ops.OnACID = 1;
    oa_tag = ''; % _OA
else
    ops.OnACID = 0;
    oa_tag = '';
end

ops.num_volt_in_files = length(ops.files_volt_in);

ops.file_path_full_ca = cell(ops.num_planes,1); 
if strcmp(ops.ca_processing ,'onacid')
    if ops.num_planes > 1
        for n_pl = 1:ops.num_planes
            ops.file_path_full_ca{n_pl} = [ops.file_dir '\' ops.file_core oa_tag '_mpl' num2str(n_pl)];
        end
    else
        ops.file_path_full_ca{1} = [ops.file_dir '\' ops.file_core oa_tag];
    end
else
    ops.file_path_full_ca = {[ops.file_dir '\' 'TracesRAW_' ops.file_core '_REDCHAN.csv']};
end

% ops.file_path_full_voltxml = cell(ops.num_volt_in_files,1);
% for n_file = 1:ops.num_volt_in_files
%     ops.file_path_full_voltxml{n_file} = [ops.file_dir '\' ops.files_volt_in{n_file}];
% end

ops.file_save_path_full = [ops.file_dir '\' ops.file_core];
% if ops.OnACID
%     ops.file_save_path_full = [ops.file_save_path_full '_OA'];
% end

ops.file_save_path_full_processing_params = [ops.file_save_path_full '_S1_processing_params.mat'];

% collect errors in this
ops.errors = {};

ops.volt_chan_labels = {'stim type', 'LED', 'Locomotion', 'TDT audio volt'};

% identify rest traces
ops.is_rest = zeros(ops.num_volt_in_files,1);
for n_file = 1:ops.num_volt_in_files
    if ~isempty(strfind(ops.files_volt_in{n_file},'rest'))
        ops.is_rest(n_file) = 1;
    else
        ops.is_rest(n_file) = 0;
    end
end

end