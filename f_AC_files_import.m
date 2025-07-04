function ops = f_AC_files_import(ops)

%% load data speciefied in xlsx file from given directory 

% drive4 = 'F';
% drive5 = 'G';
% grive6 = 'E';

% ops.file_dir = 'E:\data\AC\AC_data_OA_3_16_20';
% AC_data = readtable('AC_data_list.xlsx');
% ops.paradigm_type = 'ammn'; 

%ops.file_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\AC_data_OA_3_16_20';

idx1 = strcmpi(ops.experiment_type, {ops.experiments.name});
ops.file_dir = ops.experiments(idx1).data_path;
fprintf('Loading %s data from %s...\n', ops.experiment_type, ops.file_dir);

fpath = [ops.data_dir, ops.data_list_fname];
vars_list = {'DOB', 'Cran_date'};
opts = detectImportOptions(fpath);
for n_var = 1:numel(vars_list)
    opts = setvaropts(opts,vars_list{n_var},'InputFormat','MM/dd/uuuu'); 
end

AC_data = readtable([ops.data_dir, ops.data_list_fname], opts);

ops.paradigm_type = ops.experiments(idx1).paradigm; % 'ammn' 'freq_grating' 'cont'

ops.fpath_reg_data = [ops.gui_save_dir '\' ops.experiments(idx1).save_reg_fname];
ops.fpath_mat_data = [ops.gui_save_dir '\' ops.experiments(idx1).save_mat_fname];

%%
if numel(ops.paradigm_type)
    temp_idx = strcmpi(AC_data.paradigm, ops.paradigm_type);
    AC_data = AC_data(temp_idx,:);
end

if numel(ops.experiment_type)
    temp_idx = strcmpi(AC_data.experiment, ops.experiments(idx1).experiment);
    AC_data = AC_data(temp_idx,:);
end

AC_data = AC_data(AC_data.use_dset ~= 0,:);

%% 

num_dsets = size(AC_data,1);
f_names = cell(num_dsets,1);
f_names_dir = cell(num_dsets,1);
for n_dset = 1:num_dsets
    f_names{n_dset} = sprintf('%s_im%d_%s_%s', AC_data.mouse_id{n_dset}, AC_data.im_num(n_dset), AC_data.dset_name{n_dset}, AC_data.mouse_tag{n_dset});
    dir_data = dir([ops.file_dir '\*' AC_data.dset_name{n_dset} '*' AC_data.mouse_tag{n_dset} '*_sort.mat']);
    f_names_dir{n_dset} = {dir_data.name};
end

AC_data.dset_name_full = f_names;
AC_data.dset_name_load = f_names_dir;

%%
% use_dset = AC_data.im_use_dset;
% use_dset(isnan(use_dset)) = 0;
% 
% AC_data.mpl(isnan(AC_data.mpl)) = 0;
% 
% %use_dset(AC_data.mpl<2) = 0;
% %use_dset(AC_data.mpl>1) = 0;
% 
% AC_data = AC_data(logical(use_dset),:);
% AC_data = AC_data(strcmpi(AC_data.paradigm,ops.paradigm_type),:);

ops.conditions = unique(AC_data.area, 'stable');
ops.regions_to_analyze = unique(AC_data.area, 'stable');

% ops.file_names.A1 = AC_data.experiment(strcmpi(AC_data.area, 'A1'));
% ops.file_names.AAF = AC_data.experiment(strcmpi(AC_data.area, 'AAF'));
% ops.file_names.A2 = AC_data.experiment(strcmpi(AC_data.area, 'A2'));
% ops.file_names.UF = AC_data.experiment(strcmpi(AC_data.area, 'UF'));

ops.AC_data = AC_data;

end