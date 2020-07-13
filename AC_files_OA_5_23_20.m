
ops.file_dir = 'C:\Users\ys2605\Desktop\AC_data\AC_data_OA_3_16_20';

if ops.blah == 1
    ops.paradigm_type = 'ammn';
    
    % caiman OA datasets
     files_ammn.A1 = {'1_10_2_18',  '2_10_2_18',...
                      '1_11_21_18', '2_11_21_18',...
                      '1_5_14_20a', '2_5_14_20a',...
                      '1_5_21_20',  '2_5_21_20'...
                      };
     %, '3_10_2_18_OA', '4_10_2_18_OA', '1_10_17_18_OA', '2_10_17_18_OA',
     %'3_11_21_18_OA', '4_11_21_18_OA'};% 
     files_ammn.A2 = {'1_10_17_18', '2_10_17_18', '3_10_17_18', '4_10_17_18',...
                      '3_5_15_20',...
                      '1_5_17_20',  '2_5_17_20',...
                      '1_5_21_20',  '2_5_21_20'...
                      }; % 
     files_ammn.AAF = {'1_10_2_18', '2_10_2_18'};
     files_ammn.DF = {'1_10_16_18', '2_10_16_18', '3_10_16_18', '4_10_16_18'};
     
    ops.file_names = files_ammn;
    clear files_ammn;
elseif ops.blah == 2
    ops.paradigm_type = 'freq_grating';
    
    % caiman OA datasets
    files_freq_grating.A1 = {'1_10_2_18_OA', '2_10_2_18_OA', '3_10_2_18_OA', '4_10_2_18_OA'};
    files_freq_grating.A2 = {'1_10_17_18_OA'};
    files_freq_grating.AAF = {};
    files_freq_grating.DF = {'1_10_16_18_OA', '2_10_16_18_OA', '3_10_16_18_OA', '4_10_16_18_OA'};
    
    ops.file_names = files_freq_grating;
    clear files_freq_grating;
    
end