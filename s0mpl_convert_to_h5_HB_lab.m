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
params.data_dir = 'C:\Users\ys2605\Desktop\';
%data_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_24_21_pt3\';

% % type 2 and 3
% load_file_name = 'rest1_5_9_19.hdf5'; % only for 2 and 3
params.fname =  'some_vid2.h5'; %

%save_dir = 'C:\Users\rylab_dataPC\Desktop\Yuriy\DD_data\proc_data';
%save_dir = 'E:\data\Auditory\caiman_out_multiplane';
%save_dir = 'J:\mouse\backup\2018\caiman_out_dLGN';
%save_dir = 'L:\data\Auditory\caiman_out';
params.save_dir = 'C:\Users\ys2605\Desktop\';

params.im_target_fname = '';%'A1_cont_0.5_12_4_21a_h5cutsdata.mat';

params.num_planes = 1; % number of planes
params.do_moco = 1;
params.do_bidi = 0;

params.use_prairie_mpl_tags = 0;

params.load_type = 3; 
% 1 = Prairie tiffs
% 2 = tiff stack (needs file name)
% 3 = h5 stack (needs file name)

%%
f_s0mpl_convert_to_h5_HB_lab(params);

%%