clear;
close all; 
addpath([pwd '\s1_functions']);

% %%
% data_dir = {'H:\data\AC\2021\',...
%             'H:\data\AC\2022\'};
%          
data_dir = {'D:\data\AC\2p\2020';...
            'I:\mouse\auditory\2018'};

%data_dir = {'G:\data\Auditory\2018'};

% data_dir = {'F:\AC_data\'};
        
save_dir = {'F:\AC_data\caiman_data_missmatch\'};%,...
%save_dir = {'F:\AC_data\caiman_data_dream\'};

params.dset_table_fpath = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\AC_data_list_all.xlsx';

params.limit.dset_name =        '';
params.limit.experiment =       'missmatch';
params.limit.mouse_id =         'M10';
params.limit.mouse_tag =        '';
params.limit.dset_name =        '';
params.limit.FOV_num =          1;

%%
AC_data = f_s0_parse_tab_data(params);

mouse_id_all = unique(AC_data.mouse_id, 'stable');

%% set default params
AC_data.do_mc(isnan(AC_data.do_bidi)) = 1;
AC_data.do_bidi(isnan(AC_data.do_bidi)) = 0;
AC_data.mc_zero_edge(isnan(AC_data.mc_zero_edge)) = 1;
AC_data.mc_rigid_met(isnan(AC_data.mc_rigid_met)) = 1;
AC_data.mc_do_nonrigid(isnan(AC_data.mc_do_nonrigid)) = 0;
AC_data.mc_nonrigid_met(isnan(AC_data.mc_nonrigid_met)) = 1;
%%

if iscell(save_dir)
    params.save_dir = save_dir{1};
else
    params.save_dir = save_dir;
end

fprintf('Running %d dsets total...\n', size(AC_data,1))

for n_ms = 1:numel(mouse_id_all)
    AC_data2 = AC_data(strcmpi(AC_data.mouse_id, mouse_id_all{n_ms}),:);
    fprintf('Mouse id %s; %d dsets...\n', mouse_id_all{n_ms}, size(AC_data2,1));
    % check if folder exists
    
    % set moco target as first in list
    idx2 = logical(sum(AC_data2.im_num == unique(AC_data2.mc_to_dset)',2));
    AC_data2 = [AC_data2(idx2,:); AC_data2(~idx2,:)];
    
    for n_dset = 1:size(AC_data2,1)
        do_s0 = true;
        
        fold_name = sprintf('%s_%s_%s', AC_data2.mouse_id{n_dset}, AC_data2.mouse_tag{n_dset}, AC_data2.experiment{n_dset});
        
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
            cdset = AC_data2(n_dset,:);
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
                params.do_moco = cdset.do_mc;
                params.do_nonrigid = cdset.mc_do_nonrigid;
                params.moco_zero_edge = cdset.mc_zero_edge;
                params.do_bidi = cdset.do_bidi;
                params.moco_rigid_method = cdset.mc_rigid_met;
                params.moco_nonrigid_method = cdset.mc_nonrigid_met;
                params.dset_name = cdset.dset_name{1};
                params.save_dir;
                if or(cdset.align_pulse_crop_method == 0, cdset.align_pulse_crop_method == 2)
                    params.align_pulse_crop_method = cdset.align_pulse_crop_method;
                else
                    params.align_pulse_crop_method = 1; % default is 1 = auto; 2 = manual; 0 = full movie
                end
                
                params.im_target_fname = '';
                if params.do_moco
                    if ~isempty(cdset.mc_to_dset)
                        if cdset.im_num ~= cdset.mc_to_dset
                            source_dset = cdset.mc_to_dset;
                            source_dset_idx = AC_data2.im_num == source_dset;
                            fname_dset1 = sprintf('%s_im%d_%s_%s', AC_data2.mouse_id{source_dset_idx}, AC_data2.im_num(source_dset_idx), AC_data2.dset_name{source_dset_idx}, AC_data2.mouse_tag{source_dset_idx});
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