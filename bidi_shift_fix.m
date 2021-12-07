folder_path = 'F:\training\vincent\7_16_19\plane1_total_without_pulse.hdf5';

if ~exist('data', 'var')
    disp('Loading movie...')
    data = bigread3(folder_path, 1);
end


data = h5read(folder_path,'/mov');


[Ly, Lx, NT] = size(data);
n_frames = NT;


[data2, bidi_phase2] = f_bidi_shift(data, NT, 1);




[path, name, ~] = fileparts(folder_path);
save_file_name = [path '\' name '_shift.hdf5'];
f_save_mov_YS(data2, save_file_name, '/mov')



figure;
imagesc(sum(data(:,:,1:NT),3));
title('Pre');

figure;
imagesc(sum(data2(:,:,1:NT),3));
title('Post');


figure;
plot(bidi_phase2);

