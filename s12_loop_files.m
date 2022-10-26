clear;
close all;

%%
addpath([pwd '\s1_functions']);
addpath([pwd '\general_functions'])

%%
%ops.file_dir = 'F:\AC_data\caiman_data_dream\preprocessing';
%ops.file_dir = 'F:\AC_data\caiman_data_missmatch\preprocessing';
ops.file_dir = 'F:\AC_data\caiman_data_echo\preprocessing';

params.dset_table_fpath = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\AC_data_list_all.xlsx';

params.limit.dset_name =        '';
params.limit.experiment =       'echo';
params.limit.mouse_id =         '';
params.limit.mouse_tag =        '';
params.limit.dset_name =        '';
params.limit.FOV_num =          '';
params.limit.paradigm =         '';

%%
AC_data = f_s0_parse_tab_data(params);

mouse_id_all = unique(AC_data.mouse_id, 'stable');

%%
fprintf('Running %d dsets total...\n', size(AC_data,1))

for n_ms = 1:numel(mouse_id_all)
    AC_data2 = AC_data(strcmpi(AC_data.mouse_id, mouse_id_all{n_ms}),:);
    fprintf('Mouse id %s; %d dsets...\n', mouse_id_all{n_ms}, size(AC_data2,1));
    % check if folder exists
    
    for n_dset = 1:size(AC_data2,1)
        
        cdata = AC_data2(n_dset,:);
        ops1 = ops;
        ops1.file_core = sprintf('%s_im%d_%s_%s', cdata.mouse_id{1}, cdata.im_num, cdata.dset_name{1}, cdata.mouse_tag{1});
        
        ops1.files_volt_in = {['' ops1.file_core '_prairie']};
        
        do_s12 = 1;
        if exist([ops1.file_dir '\' ops1.file_core '_processed_data.mat'], 'file')
            do_s12 = 0;
            fprintf('skipping, already exists: ');
        else
            fprintf('Running: ');
        end
        fprintf('%s\n', ops1.file_core);
        
        if ~isnan(cdata.s12_volt_chan)
            ops1.volt_chan_order = str2double(num2cell(num2str(cdata.s12_volt_chan)));
        else
            ops1.volt_chan_order = [1 2 3 4];
        end
        
        if ~isnan(cdata.s12_volt_chan_bh)
            ops1.volt_chan_order_bh = str2double(num2cell(num2str(cdata.s12_volt_chan_bh)));
        end
        
        if ~isnan(cdata.s12_volt_chan_stim)
            ops1.volt_chan_order_stim = str2double(num2cell(num2str(cdata.s12_volt_chan_stim)));
        end
        
        if ~isnan(cdata.s12_exp_win_sel)
            ops1.exp_win_selection = cdata.s12_exp_win_sel;
        end
        
        ops1.paradigm = cdata.paradigm{1};
        
        ops1.num_planes = cdata.mpl;
        
        if do_s12
            f_s12_preprocess_voltage_ca_data(ops1);
        end
    end
end

fprintf('All done\n')