
mov_dir = 'F:\AC_data\caiman_data_echo\movies\';
%mov_dir = 'G:\data\Auditory\caiman_data_missmatch\movies';

fname = {'M226_im9_A2_cont44_7_17_22_mpl5_pl1.h5',...
         %'M10_im12_A1_ammn2_5_31_20.h5',...
         %'M10_im13_A1_freq_grating1_5_31_20.h5',...
         %'M10_im14_A1_freq_grating2_5_31_20.h5',...
         };

 max_size = 25000;
 
for n_fl = 1:numel(fname)
    [~, fname2, ext1] = fileparts(fname{n_fl});

    Y = h5read([mov_dir '\' fname2 ext1],'/mov');

    [d1, d2, T] = size(Y);
    
    num_blocks = ceil(T/max_size);
    start1 = 1;
    for n_bl = 1:num_blocks
        end1 = min((start1 + max_size - 1), T);
        f_save_mov_YS(Y(:,:,start1:end1), sprintf('%s\\%s_pt%d%s', mov_dir, fname2, n_bl, ext1), '/mov');
        start1 = end1 + 1;
    end
end

