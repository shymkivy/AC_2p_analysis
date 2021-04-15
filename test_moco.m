clear;
close all;


pwd1 = fileparts(which('test_moco.m'));
addpath([pwd1 '\general_functions']);

mov_load_path = 'F:\data\Auditory\2018\10_2_18_im\A2_freq_grating-001';

mpl_tags = {'Ch2'};

save_path = 'F:\data\Auditory\caiman_out\movies\';

Y = f_collect_prairie_tiffs4(mov_load_path, mpl_tags{1});


Y_reg = f_register_suite2p_YS(Y);

f_save_mov_YS(Y, sprintf('%sA2_freq_grating1_s2p_MC_pre.h5', save_path), '/mov')
f_save_mov_YS(Y_reg, sprintf('%sA2_freq_grating1_s2p_MC-post.h5', save_path), '/mov')