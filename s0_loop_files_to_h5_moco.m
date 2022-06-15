clear;
close all;
addpath([pwd '\s1_functions']);

% %%
%  data_dir = {'D:\data\AC\2022\',...
%              'D:\data\AC\2022\'};
         
% data_dir = {'G:\data\Auditory\2018\',...
%             'E:\data\AC\2p\2020\'};
 
data_dir = {'F:\AC_data\'};
        
%save_dir = {'F:\AC_data\caiman_data_missmatch\'};%,...
save_dir = {'F:\AC_data\caiman_data_dream3\'};

params.dset_table_fpath = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\AC_data_list_all.xlsx';

experiment_tag = 'dream';

limilt_mouse_id = 'M125';
limit_mouse_tag = '';
limit_dset_name = 'AC_ammn4';

%%
AC_data = readtable(params.dset_table_fpath);

%%

AC_data = AC_data(~isnan(AC_data.idx),:);

AC_data = AC_data(AC_data.do_proc == 1,:);

%%
AC_data2 = AC_data(strcmpi(AC_data.experiment, experiment_tag),:);

AC_data3 = AC_data2;
if exist('limilt_mouse_id', 'var')
    if numel(limilt_mouse_id)
        AC_data3 = AC_data3(strcmpi(AC_data3.mouse_id, limilt_mouse_id),:);
    end
end

if exist('limit_mouse_tag', 'var')
    if numel(limit_mouse_tag)
        AC_data3 = AC_data3(strcmpi(AC_data3.mouse_tag, limit_mouse_tag),:);
    end
end

if exist('limit_dset_name', 'var')
    if numel(limit_dset_name)
        AC_data3 = AC_data3(strcmpi(AC_data3.dset_name, limit_dset_name),:);
    end
end

mouse_id_all = unique(AC_data3.mouse_id, 'stable');

%% set default params
AC_data3.do_moco(isnan(AC_data3.do_bidi)) = 1;
AC_data3.do_bidi(isnan(AC_data3.do_bidi)) = 0;

AC_data3.moco_zero_edge(isnan(AC_data3.moco_zero_edge)) = 1;

%%

if iscell(save_dir)
    params.save_dir = save_dir{1};
else
    params.save_dir = save_dir;
end

fprintf('Running %d dsets total...\n', size(AC_data3,1))

for n_ms = 1:numel(mouse_id_all)
    AC_data4 = AC_data3(strcmpi(AC_data3.mouse_id, mouse_id_all{n_ms}),:);
    fprintf('Mouse id %s; %d dsets...\n', mouse_id_all{n_ms}, size(AC_data4,1));
    % check if folder exists
    
    for n_dset = 1:size(AC_data4,1)
        do_s0 = true;
        
        fold_name = sprintf('%s_%s_%s', AC_data4.mouse_id{n_dset}, AC_data4.mouse_tag{n_dset}, experiment_tag);
        
        if isstring(data_dir)
            data_dir = {data_dir};
        end
        
        fold_exist = false(numel(data_dir),1);
        for n_data_dir = 1:numel(data_dir)
            data_dir2 = data_dir{n_data_dir};
            if exist([data_dir2 '\' fold_name], 'dir')
                fold_exist(n_data_dir) = 1;
            end
        end
        
        if ~sum(fold_exist)
            do_s0 = 0;
            warning(['Data directory does not exist: ' fold_name])
        else
            params.data_dir = [data_dir{fold_exist} '\' fold_name];
        end
        
        if do_s0
            cdset = AC_data4(n_dset,:);
            params.fname = sprintf('%s_im%d_%s_%s', cdset.mouse_id{1}, cdset.im_num, cdset.dset_name{1}, cdset.mouse_tag{1});

            % check it output already exists
            num_match = 0;
            if iscell(save_dir)
                for n_dir = 1:numel(save_dir)
                    dir_list = dir([save_dir{n_dir} '\movies\*.h5']);
                    dir_names = {dir_list.name};
                    for n_file = 1:numel(dir_names)
                        pat1 = strfind(dir_names{n_file}, params.fname);
                        if ~isempty(pat1)
                            num_match = num_match + 1;
                        end
                    end
                end
            else
                dir_list = dir([save_dir{n_dir} '\movies\*.h5']);
                dir_names = {dir_list.name};
                for n_file = 1:numel(dir_names)
                    pat1 = strfind(dir_names{n_file}, params.fname);
                    if ~isempty(pat1)
                        num_match = num_match + 1;
                    end
                end
            end

            if ~num_match
                params.num_planes = cdset.mpl;
                params.do_moco = cdset.do_moco;
                params.moco_zero_edge = cdset.moco_zero_edge;
                params.do_bidi = cdset.do_bidi;
                params.dset_name = cdset.dset_name{1};
                params.save_dir;
                if or(cdset.align_pulse_crop_method == 0, cdset.align_pulse_crop_method == 2)
                    params.align_pulse_crop_method = cdset.align_pulse_crop_method;
                else
                    params.align_pulse_crop_method = 1; % default is 1 = auto; 2 = manual; 0 = full movie
                end
                
                params.im_target_fname = '';
                if params.do_moco
                    if ~isempty(cdset.moco_to_dset)
                        if cdset.im_num ~= cdset.moco_to_dset
                            source_dset = cdset.moco_to_dset;
                            source_dset_idx = AC_data4.im_num == source_dset;
                            fname_dset1 = sprintf('%s_im%d_%s_%s', AC_data4.mouse_id{source_dset_idx}, AC_data4.im_num(source_dset_idx), AC_data4.dset_name{source_dset_idx}, AC_data4.mouse_tag{source_dset_idx});
                            params.im_target_fname = fname_dset1;
                            fprintf('Using target registration im from %s\n', params.im_target_fname);
                        end
                    end
                end
                f_s0mpl_convert_to_h5(params);
            else
                fprintf('%s already exists, moving on...\n', params.fname)
            end

        end
    end
end

fprintf('All done\n')