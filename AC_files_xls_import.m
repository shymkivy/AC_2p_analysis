%% load data speciefied in xlsx file from given directory 

% drive4 = 'F';
% drive5 = 'G';
% grive6 = 'E';

% ops.file_dir = 'E:\data\AC\AC_data_OA_3_16_20';
% AC_data = readtable('AC_data_list.xlsx');
% ops.paradigm_type = 'ammn'; 

ops.file_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';
AC_data = readtable('AC_data_list_echo.xlsx');
ops.paradigm_type = 'cont'; % 'ammn' 'freq_grating' 'cont'

%%
use_dset = AC_data.im_use_dset;
use_dset(isnan(use_dset)) = 0;

AC_data.mpl(isnan(AC_data.mpl)) = 0;

%use_dset(AC_data.mpl<2) = 0;
%use_dset(AC_data.mpl>1) = 0;

AC_data = AC_data(logical(use_dset),:);
AC_data = AC_data(strcmpi(AC_data.paradigm,ops.paradigm_type),:);

ops.file_names.A1 = AC_data.experiment(strcmpi(AC_data.area, 'A1'));
ops.file_names.AAF = AC_data.experiment(strcmpi(AC_data.area, 'AAF'));
ops.file_names.A2 = AC_data.experiment(strcmpi(AC_data.area, 'A2'));
ops.file_names.UF = AC_data.experiment(strcmpi(AC_data.area, 'UF'));

ops.AC_data = AC_data;
clear AC_data use_dset;

