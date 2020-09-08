
ops.file_dir = 'E:\data\AC\AC_data_OA_3_16_20';

if ops.blah == 1
    ops.paradigm_type = 'ammn';
    
    % caiman OA datasets
     files_ammn.A1 = {'1_10_2_18',  '2_10_2_18',...  %'3_10_2_18', '4_10_2_18', 
                      '1_10_17_18', '2_10_17_18',... %'3_10_17_18', '4_10_17_18',
                      '1_11_21_18', '2_11_21_18',... %'3_11_21_18', '4_11_21_18'
                      '1_5_14_20a', '2_5_14_20a',...
                      '1_5_21_20',  '2_5_21_20',...
                      '3_5_25_20',  '4_5_25_20',...
                      '1_5_31_20',  '2_5_31_20',...
                      };
     files_ammn.A2 = {'1_10_17_18', '2_10_17_18', '3_10_17_18', '4_10_17_18',...
                      '3_5_15_20',...   % blurry
                      '1_5_17_20',  '2_5_17_20',...
                      '1_5_21_20',  '2_5_21_20',...
                      '1_5_25_20',  '2_5_25_20',...
                      '1_5_31_20',  '2_5_31_20',...
                      };
     % to OA: '2_10_2_18',
     % to OA: '4_5_15_20'
     files_ammn.AAF = {'1_10_2_18', '2_10_2_18',... % '3_10_2_18', '4_10_2_18'
                       '1_5_15_20', '2_5_15_20',...
                       '1_5_17_20', '2_5_17_20',...
                       '1_5_21_20', '2_5_21_20',...
                       '1_5_25_20', '2_5_25_20',...
                       '1_5_31_20', '2_5_31_20', ...
                      };
                  
     files_ammn.DF = {'1_10_16_18', '2_10_16_18',... '3_10_16_18', '4_10_16_18'...
                      '1_5_15_20', '2_5_15_20',...
                      '1_5_21_20', '2_5_21_20',...
                      '1_5_25_20', '2_5_25_20',...
                      '1_5_31_20', '2_5_31_20',...
                     };
     % to OA: '1_5_14_20a', '2_5_14_20a'
    ops.file_names = files_ammn;
    clear files_ammn;
elseif ops.blah == 2
    ops.paradigm_type = 'freq_grating';
    
    % caiman OA datasets
    files_freq_grating.A1 = {'1_10_2_18',  '2_10_2_18',... %'3_10_2_18', '4_10_2_18',...
                             '1_10_17_18', '2_10_17_18',...
                             '1_11_21_18', '2_11_21_18',... %'3_11_21_18',  '4_11_21_18',...
                             '1_5_14_20a', '2_5_14_20a',...
                             '1_5_21_20', ...
                             '3_5_25_20',  '4_5_25_20',...
                             '1_5_31_20', ...
                             };
        % to do: '2_5_21_20',  '2_5_31_20'
                         
    files_freq_grating.A2 = {'1_10_17_18', '2_10_17_18',...
                             '1_5_17_20', '2_5_17_20',...
                             '1_5_21_20',...
                             '1_5_25_20', '2_5_25_20',...
                             '1_5_31_20', '2_5_31_20',...
                             };
        % TO ADD : '2_5_21_20', '3_5_15_20', '4_5_15_20', '3_10_17_18', '4_10_17_18'
    files_freq_grating.AAF = {'1_5_15_20', '2_5_15_20',...
                              ...
                              '1_5_21_20',...
                              '2_5_25_20',...
                              '1_5_31_20', '2_5_31_20',...
                              };
        % to do: '1_5_17_20','2_5_17_20', '2_5_17_20', '2_5_21_20', '1_5_25_20', 
    files_freq_grating.DF = {'1_10_16_18', '2_10_16_18', '3_10_16_18', '4_10_16_18'...
                              ...
                             '1_5_15_20',  '2_5_15_20',...
                             '1_5_21_20',  '2_5_21_20',...
                             '1_5_25_20',  '2_5_25_20',...
                             '1_5_31_20',  '2_5_31_20',...
                             };
        % to do: '1_5_14_20a', '2_5_14_20a',
    ops.file_names = files_freq_grating;
    clear files_freq_grating;
    
end