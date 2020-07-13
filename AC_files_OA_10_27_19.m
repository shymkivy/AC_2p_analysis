
ops.file_dir = 'C:\Users\ys2605\Desktop\A1_data\AC_data_OA_10_27_19';

if ops.blah == 1
    ops.paradigm_type = 'ammn';
    
    % caiman OA datasets
     files_ammn.A1 = {'1_10_2_18_OA', '2_10_2_18_OA', '3_10_2_18_OA', '4_10_2_18_OA', '1_10_17_18_OA', '2_10_17_18_OA'};
     files_ammn.A2 = {};
     files_ammn.AAF = {'1_10_2_18_OA', '2_10_2_18_OA', '3_10_2_18_OA', '4_10_2_18_OA'};
     files_ammn.DAM = {'1_10_16_18_OA', '2_10_16_18_OA', '3_10_16_18_OA', '4_10_16_18_OA'};
     
    ops.file_names = files_ammn;
    clear files_ammn;
elseif ops.blah == 2
    ops.paradigm_type = 'freq_grating';
    
    % caiman OA datasets
    files_freq_grating.A1 = {'1_10_2_18_OA', '2_10_2_18_OA', '3_10_2_18_OA', '4_10_2_18_OA'};
    files_freq_grating.A2 = {};
    files_freq_grating.AAF = {};
    files_freq_grating.DAM = {'1_10_16_18_OA', '2_10_16_18_OA', '3_10_16_18_OA', '4_10_16_18_OA'};
    
    ops.file_names = files_freq_grating;
    clear files_freq_grating;
    
end