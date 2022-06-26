clear;
close all;

addpath('C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\s1_functions');
addpath('C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\general_functions');

data_dir = 'D:\data\caiman_data_dream';
%data_dir = 'F:\AC_data\caiman_data_dream3';
dset_name = 'M125_im7_AC_ammn_stim7_4_26_22';

save_dir = 'D:\data\AC\2021\M125_4_26_22_dream';
%save_dir = 'F:\AC_data\M108_2_4_22a_dream';


%stim_chan = 7;

%%

dir_mov = [data_dir '\movies\'];
dir_proc = [data_dir '\preprocessing\'];

%%
f_list = dir([dir_mov '*' dset_name '*.h5']);

num_planes = numel(f_list);

Y = cell(num_planes,1);
for n_pl = 1:num_planes
    Y{n_pl} = h5read([dir_mov f_list(n_pl).name], '/mov');
end

%% load proc data

proc_data = load([dir_proc dset_name '_processed_data.mat']);

if exist([dir_proc dset_name '_mpl_scan.mat'], 'file')
    data_mplscan = load([dir_proc dset_name '_mpl_scan.mat']);
    if isfield(data_mplscan.scan_data, 'custom_stim_data')
        same_stim = 1;
    else
        same_stim = 0;
    end
else
    same_stim = 0;
end

%data_h5cuts = load([dir_proc dset_name '_h5cutsdata.mat']);

%data_volt = dlmread([dir_proc dset_name '_prairie.csv'],',',1,0);

%data_xml = extract_frame_data_from_XML2([dir_proc dset_name '_prairie.xml']);

%% fill vid cuts

[d1, d2, num_frames_cut] = size(Y{1});

vid_cuts = proc_data.data.file_cuts_params{1}.vid_cuts_trace;

num_frames = numel(vid_cuts);

if num_frames > num_frames_cut
    for n_pl = 1:num_planes
        Y_temp = zeros(d1, d2, num_frames, 'uint16');
        Y_temp(:,:,logical(vid_cuts)) = Y{n_pl};
        Y{n_pl} = Y_temp;
    end
end

%%

onset_times_plockle = proc_data.data.stim_times_volt{strcmpi(proc_data.ops.chan_labels, 'Pockel')};
%onset_times_stim_trig = volt_time_s(f_get_pulse_times(volt_stim_trig_n,0.5, 99999) - 1);
frame_times = proc_data.data.frame_data.frame_times_mpl{1};

if same_stim
    stim_idx = data_mplscan.scan_data.custom_stim_data.stim_patterns_index;
    stim_pat_id = data_mplscan.scan_data.custom_stim_data.stim_patterns_id;
    stim_tab = data_mplscan.scan_data.group_table_stim;
else
    num_patterns = 1;
end


stim_resp_frames = 20;
base_frames = 10;

for n_pat = 1 % 31:50
    if same_stim
        pat_tab = stim_tab(stim_tab.Pattern == stim_pat_id(n_pat),:);
        stim_times = onset_times_plockle(stim_idx == n_pat);
    else
        stim_times = onset_times_plockle;
    end
    stim_times(stim_times<base_frames) = [];
    
    
    num_stim = numel(stim_times);

    stim_times_frame = zeros(num_stim,1);
    for n_st = 1:num_stim
        stim_times_frame(n_st) = (find(stim_times(n_st) < frame_times,1) + 1);
    end

    frames_base_im = cell(num_planes, num_stim);
    frames_im = cell(num_planes, num_stim);
    frames_all = cell(num_planes, num_stim);
    for n_st = 1:num_stim
        for n_pl = 1:num_planes
            frames_im{n_pl, n_st} = mean(Y{n_pl}(:,:,stim_times_frame(n_st):(stim_times_frame(n_st)+stim_resp_frames)),3);
            frames_base_im{n_pl, n_st} = mean(Y{n_pl}(:,:,(stim_times_frame(n_st)-base_frames):(stim_times_frame(n_st)-1)),3);
            frames_all{n_pl, n_st} = Y{n_pl}(:,:,(stim_times_frame(n_st)-20):(stim_times_frame(n_st)+30));
            
        end
    end


    frame_st_ave = cell(num_planes,1);
    frame_st_ave_base = cell(num_planes,1);
    framse_all_ave_cat = cell(num_planes,1);
    for n_pl = 1:num_planes
        frame_st_ave{n_pl} = mean(cat(3,frames_im{n_pl,:}),3);
        frame_st_ave_base{n_pl} = mean(cat(3,frames_base_im{n_pl,:}),3);
        framse_all_ave_cat{n_pl} = mean(cat(4,frames_all{n_pl,:}),4);
    end
    frame_st_ave_mpl = cat(2,frame_st_ave{:});
    frame_st_ave_mpl_base = cat(2,frame_st_ave_base{:});
    framse_all_all_cat_mpl = cat(2,framse_all_ave_cat{:});
    
    f_save_tif_stack2_YS(framse_all_all_cat_mpl, sprintf('%s\\stim_resp_pat%d.tiff', save_dir, n_pat));
    
    figure; 
    subplot(3,1,1); axis equal tight;
    imagesc(frame_st_ave_mpl);
    subplot(3,1,2); axis equal tight;
    imagesc(frame_st_ave_mpl_base);
    subplot(3,1,3); axis equal tight;
    imagesc(frame_st_ave_mpl - frame_st_ave_mpl_base);
    %sgtitle(sprintf('cell stim %d; x=%.1f, y=%.1f, z=%d', n_cell, pat_tab.X(1), pat_tab.Y(1), pat_tab.Z(1)));
end


for n_st = 1:10
    frame_full = cat(2, frames_im{:,n_st});
    figure; imagesc(frame_full); axis equal tight;
    title(sprintf('stim %d', n_st))
end

figure; imagesc(mean(frames_im,3))



%%
zoom = 1.2;
foc_size = 511;
z_range = 1;

gr_data_stim = data_mplscan.scan_data.group_table_stim;

pl_ave = cell(num_planes,1);
for n_pl = 1:num_planes
    pl_ave{n_pl} = mean(Y{n_pl},3);
end

pl_ave_all = cat(2, pl_ave{:});

figure;
imagesc(pl_ave_all);
axis equal tight;
title(sprintf('planes 1 - %d', num_planes))
%%

figure; imagesc(mean(Y{1},3))
figure; imagesc(mean(Y{2},3))




