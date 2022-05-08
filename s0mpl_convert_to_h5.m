%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Workflow
%       Load movie, 3 options
%           1: Prairie tiffs
%           2: Tiff stack
%           3: H5 stack
%       Crop pulses
%       Save as H5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear;
close all;
addpath([ pwd '\s1_functions']);

%%
params.data_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\M125_bh_5_2_22\';
%data_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_24_21_pt3\';

save_prefix = 'M125_im10_';
fname = 'AC_tone_lick_reward'; %
file_num = '1';
file_date = '5_2_22';
%file_date = '10_2_18';
% % type 2 and 3
% load_file_name = 'rest1_5_9_19.hdf5'; % only for 2 and 3
params.fname = [save_prefix, fname, file_num '_' file_date];
params.dset_name = [fname, file_num];

%save_dir = 'C:\Users\rylab_dataPC\Desktop\Yuriy\DD_data\proc_data';
%save_dir = 'E:\data\Auditory\caiman_out_multiplane';
%save_dir = 'J:\mouse\backup\2018\caiman_out_dLGN';
%save_dir = 'L:\data\Auditory\caiman_out';
params.save_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data_dream';

params.im_target_fname = '';%'A1_cont_0.5_12_4_21a_h5cutsdata.mat';

params.num_planes = 5; % number of planes or 0
params.do_moco = 1;
params.do_bidi = 0;

params.align_pulse_crop_method = 1;% default is 1 = auto; 2 = manual; 0 = full movie

%%
f_s0mpl_convert_to_h5(params);

%%