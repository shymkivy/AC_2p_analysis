ops.file_dir = 'C:\Users\ys2605\Desktop\AC_data\AC_data_OA_mpl_12_11_19';

if ops.blah == 1
    ops.paradigm_type = 'ammn';
    
    % multiplane data
    files_ammn.A1 = {};
    files_ammn.A2 = {};
    files_ammn.AAF = {};
    files_ammn.DF = {'_2_dplanes1_10_14_19', '_2_dplanes2_10_14_19'};
     
    ops.file_names = files_ammn;
    clear files_ammn;
elseif ops.blah == 2
    ops.paradigm_type = 'freq_grating';
    
    % caiman OA datasets
    files_freq_grating.A1 = {};
    files_freq_grating.A2 = {};
    files_freq_grating.AAF = {};
    files_freq_grating.DF = {};
    
    ops.file_names = files_freq_grating;
    clear files_freq_grating;
    
end