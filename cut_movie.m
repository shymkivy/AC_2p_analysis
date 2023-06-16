addpath('C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\general_functions');

fpath =  'F:\AC_data\caiman_data_missmatch\movies\';
fname = 'M8_im1_A2_ammn1_5_21_20.h5';

num_frames = 20000;

[~, fname2, ext] = fileparts(fname);

Y = h5read([fpath, fname], '/mov');

Y2 = Y(:,:,1:num_frames);

f_save_mov_YS(Y2, [fpath, fname2, '_', num2str(num_frames), ext], '/mov');


