

mov_dir = 'F:\AC_data\caiman_data_missmatch\movies\';

fname = {'M10_im1_A2_ammn1_5_31_20.h5',...
         %'M10_im2_A2_ammn2_5_31_20.h5',...
         %'M10_im3_A2_freq_grating1_5_31_20.h5',...
         %'M10_im4_A2_freq_grating2_5_31_20.h5',...
         };

for n_fl = 1:numel(fname)
    [~, fname2, ext1] = fileparts(fname{n_fl});

    Y = h5read([mov_dir '\' fname2 ext1],'/mov');

    [d1, d2, T] = size(Y);
    t05 = floor(T/2);

    f_save_mov_YS(Y(:,:,1:t05), [mov_dir '\' fname2 '_pt1' ext1], '/mov');
    f_save_mov_YS(Y(:,:,(t05+1):end), [mov_dir '\' fname2 '_pt2' ext1], '/mov');
end







