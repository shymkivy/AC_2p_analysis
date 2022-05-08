clear;
close all;

%%
addpath([pwd '\s1_functions']);
addpath([pwd '\general_functions'])

%%
ops.file_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data_dream\preprocessing';

params.dset_table_fpath = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\AC_data_list_all.xlsx';

limit_title_tag = '';

limit_experiment_tag = 'dream';
limit_mouse_id = 'M125';
limit_mouse_tag = '';

%%
AC_data = readtable(params.dset_table_fpath);

%%

AC_data = AC_data(~isnan(AC_data.idx),:);

AC_data = AC_data(AC_data.do_proc == 1,:);

if numel(limit_title_tag)
    ammns_all = false(numel(AC_data.dset_name),1);
    for n_dset = 1:numel(AC_data.dset_name)
        pat1 = strfind(AC_data.dset_name{n_dset}, limit_title_tag);
        if ~isempty(pat1)
            ammns_all(n_dset) = 1;
        end
    end
    AC_data = AC_data(ammns_all,:);
end



%%
AC_data2 = AC_data(strcmpi(AC_data.experiment, limit_experiment_tag),:);

AC_data3 = AC_data2;
if exist('limit_mouse_id', 'var')
    if numel(limit_mouse_id)
        AC_data3 = AC_data3(strcmpi(AC_data3.mouse_id, limit_mouse_id),:);
    end
end

if exist('limit_mouse_tag', 'var')
    if numel(limit_mouse_tag)
        AC_data3 = AC_data3(strcmpi(AC_data3.mouse_tag, limit_mouse_tag),:);
    end
end

AC_data3 = AC_data3(AC_data3.do_proc == 1,:);

mouse_id_all = unique(AC_data3.mouse_id, 'stable');

%%
fprintf('Running %d dsets total...\n', size(AC_data3,1))

for n_ms = 1:numel(mouse_id_all)
    AC_data4 = AC_data3(strcmpi(AC_data3.mouse_id, mouse_id_all{n_ms}),:);
    fprintf('Mouse id %s; %d dsets...\n', mouse_id_all{n_ms}, size(AC_data4,1));
    % check if folder exists
    
    for n_dset = 1:size(AC_data4,1)
        
        cdata = AC_data4(n_dset,:);
         
        ops.file_core = sprintf('%s_im%d_%s_%s', cdata.mouse_id{1}, cdata.im_num, cdata.dset_name{1}, cdata.mouse_tag{1});
        
        ops.files_volt_in = {['' ops.file_core '_prairie']};
        
        do_s12 = 1;
        if exist([ops.file_dir '\' ops.file_core '_processed_data.mat'], 'file')
            do_s12 = 0;
            fprintf('skipping: %s; already exists\n', ops.file_core)
        end
        
        if ~isnan(cdata.s12_volt_chan)
            ops.volt_chan_order = str2double(num2cell(num2str(cdata.s12_volt_chan)));
        else
            ops.volt_chan_order = [1 2 3 4];
        end
        
        ops.num_planes = cdata.mpl;
        
        if strcmpi(cdata.paradigm, 'spont') || strcmpi(cdata.paradigm, 'rest')
            ops.has_stim = 0;
        else
            ops.has_stim = 1;
        end
        
        if do_s12
            f_s12_preprocess_voltage_ca_data(ops);
        end
    end
end

fprintf('All done')