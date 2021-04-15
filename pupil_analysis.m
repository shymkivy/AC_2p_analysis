
clear;
close all;

fpath = 'C:\Users\ys2605\Desktop\stuff\AC_data\behavior_cam\3_27_21\';

%%
video_input= VideoReader([fpath, 'A2_ammn1_eye0001 21-04-11 13-13-03.avi']);

%video_input= VideoReader([fpath, 'A2_ammn1_eye_cut.mp4']);
fr1 = video_input.FrameRate;

num_frames = video_input.NumFrames;
%num_frames = 5000;
%%
frame1 = rgb2gray(video_input.readFrame);

[d1, d2] = size(frame1);

frames = zeros([d1, d2, num_frames], 'uint8');

frames(:,:,1) = frame1;

f = waitbar(0,'loading...');
for n_fr = 2:num_frames
    frames(:,:,n_fr) = rgb2gray(video_input.readFrame);
    waitbar(n_fr/num_frames,f,'loading...');
end
close(f)

%%
figure; plot(squeeze(mean(mean(frames)))); title('mean loaded trace')

%%

figure; imagesc(mean(frames,3)); axis equal tight;
title('Mean video; Select whole eye (2 clicks)')
[x_eye, y_eye] = ginput(2);
x_eye = round(x_eye);
y_eye = round(y_eye);
rectangle('Position',[min(x_eye) min(y_eye) diff(x_eye)  diff(y_eye)],'EdgeColor', 'r')
title('Mean video full')

%%
frames_cut = frames(y_eye(1):y_eye(2),x_eye(1):x_eye(2),:);

mean_frame_cut = mean(frames_cut,3);
clear frames;

%%
figure; imagesc(mean_frame_cut); hold on; axis equal tight;
title('Mean frame; Select inside pupil region for max val (2 clicks)')
[x_pup, y_pup] = ginput(2);
x_pup = round(x_pup);
y_pup = round(y_pup);
rectangle('Position',[min(x_pup) min(y_pup) diff(x_pup)  diff(y_pup)],'EdgeColor', 'g')
title('Mean frame; Now select lamp artifact (2 clicks)')
[x_lamp, y_lamp] = ginput(2);
x_lamp = round(x_lamp);
y_lamp = round(y_lamp);
rectangle('Position',[min(x_lamp) min(y_lamp) diff(x_lamp)  diff(y_lamp)],'EdgeColor', 'y')
title('Mean frame; select ceter of pupil (1 click)')
[cent_x, cent_y] = ginput(1);
cent_x = round(cent_x);
cent_y = round(cent_y);
plot(cent_x, cent_y, 'ro')
title('Mean frame')


%% play vid
% 
% num_frames2 = size(frames2,3);
% 
% figure;
% im1 = imagesc(frames2(:,:,1));
% 
% for n_fr = 2:num_frames2
%     pause(1/60);
%     im1.CData = frames2(:,:,n_fr);
% end

%%

ave_pup_int = squeeze(mean(mean(frames_cut(min(y_pup):min(y_pup)+diff(y_pup),min(x_pup):min(x_pup)+diff(x_pup),:))));

figure; plot(ave_pup_int)
title('Mean pupil intensity; select threshold for imaging onset (1 click)')
[~,thresh1] = ginput(1);
title('Mean pupil intensity')

temp_vals = find(ave_pup_int>thresh1);
onset_frame = temp_vals(1);
%% cover lamp artifact?
% for n_fr = 1:size(frames_cut,3)
%     frames_cut(min(y_lamp):(min(y_lamp)+diff(y_lamp)),min(x_lamp):(min(x_lamp)+diff(x_lamp)),n_fr) = ave_pup_int(n_fr);
% end
%%
frames_cut_df = single(frames_cut(:,:,onset_frame:end)); % 
num_frames2 = size(frames_cut_df,3);
slice_width = 3;

f = waitbar(0,'medfiltering...');
for n_fr = 1:num_frames2
    frames_cut_df(:,:,n_fr) = medfilt2(frames_cut_df(:,:,n_fr),[7, 7],'symmetric');
    waitbar(n_fr/num_frames2,f,'medfiltering...');
end
close(f)

figure; imagesc(frames_cut_df(:,:,3))

slice_y = squeeze(mean(frames_cut_df((cent_y-slice_width):(cent_y+slice_width),:,:),1));
slice_x = squeeze(mean(frames_cut_df(:,(cent_x-slice_width):(cent_x+slice_width),:),2));

slice_y_n = slice_y - min(slice_y);
slice_y_n = slice_y_n./mean(slice_y_n(x_pup(1):x_pup(2),:));

slice_x_n = slice_x - min(slice_x);
slice_x_n = slice_x_n./mean(slice_x_n(y_pup(1):y_pup(2),:));

%%
yup_loc = zeros(num_frames2,1);
ydn_loc = zeros(num_frames2,1);
xl_loc = zeros(num_frames2,1);
thresh = .7;
for n_fr = 1:num_frames2
    temp_trace = slice_x_n(cent_y:end,n_fr);
    temp_loc = find(temp_trace<thresh);
    if numel(temp_loc)
        ydn_loc(n_fr) = cent_y + temp_loc(1);
    end
    temp_trace = flipud(slice_x_n(1:cent_y,n_fr));
    temp_loc = find(temp_trace<thresh);
    if numel(temp_loc)
        yup_loc(n_fr) = cent_y - temp_loc(1);
    end
end

%%

yup_loc_f = medfilt1(yup_loc, 9);
ydn_loc_f = medfilt1(ydn_loc, 9);


%%
pupil_diam = ydn_loc - yup_loc;

pupil_diam_f = ydn_loc_f - yup_loc_f;




%%
addpath('C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\general_functions');

proc_data_fname = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data\A2_ammn1_3_27_21_processed_data.mat';
proc_data = load(proc_data_fname);

%% go from 60hz to 10hz

x = (1:num_frames2)/fr1;
xq = proc_data.data.frame_data.frame_times_mpl{1,1}/1000;

pupil_diam_ds = interp1(x,pupil_diam,xq);


figure; hold on;
plot(x, pupil_diam); axis tight;
plot(x, pupil_diam_f)
xlabel('Sec')

%%
figure;hold on;
plot(x, pupil_diam); axis tight;
plot(xq, pupil_diam_ds)
plot(xq, proc_data.data.volt_data_binned{1,1}*100)
xlabel('Sec'); ylabel('diameter (pixels)')

%%

n_pl = 3;
base_onset_win = [5 15];
%tt = 270;
tt = [170, 270, 3, 6];

col1 = {[1 .3 .3], [1 .3 .3], [.7 .7 .7], [.7 .7 .7]};
col2 = {'r', 'r', 'k', 'k'};
%%
stim_frame_index = proc_data.data.stim_frame_index{1};
trial_types = proc_data.data.trial_types;
MMN_orientations = proc_data.data.MMN_orientations;

t_p = (-4:base_onset_win(2))/10;

figure; hold on; axis tight;
pup_mean_tr = cell(4,1);
trial_means = cell(4,1);
for n_tt = 1:4
    tt_temp = tt(n_tt);
    stim_frame_index1 = stim_frame_index(trial_types == tt_temp);
    num_trials = numel(stim_frame_index1);
    
    pupil_sort = squeeze(f_get_stim_trig_resp(pupil_diam_ds', stim_frame_index1, base_onset_win));
    pupil_sort_n = pupil_sort - mean(pupil_sort(1:5,:));
    
    rem_thr = -80;
    rem_trial = false(num_trials,1);
    for n_tr = 1:num_trials
        if sum(pupil_sort_n(:,n_tr)<rem_thr)
            rem_trial(n_tr) =1;
        end
    end
    
    pupil_sort_n2 = pupil_sort_n(:,~rem_trial);
    pupil_sort2 = pupil_sort(:,~rem_trial);
    
    plot(t_p, pupil_sort2, '--', 'color', col1{n_tt})
    pup_mean_tr{n_tt} = mean(pupil_sort2,2);
    
    trial_means{n_tt} = mean(pupil_sort2);
end
for n_tt = 1:4
    plot(t_p, pup_mean_tr{n_tt}, 'Color', col2{n_tt}, 'Linewidth', 2)
end
xlabel('Time (sec)'); ylabel('Pupil Diameter (pix)')

figure; hold on;
for n_tt = 1:4
    mean2 = mean(pup_mean_tr{n_tt});
    plot(n_tt, pup_mean_tr{n_tt}, 'o', 'Color', col2{n_tt})
    
end
%% play vid with fitted data
n_fr = 1;
f1 = figure; hold on; axis equal tight; axis off
im1 = imagesc(frames_cut_df(:,:,n_fr)); set(gca,'YDir','reverse');
im_yup = plot(cent_x, yup_loc(n_fr), 'ro');
im_ydn = plot(cent_x, ydn_loc(n_fr), 'go');
pos = f1.Position;
margx = -84;
margy = -45;
rect1 = [-margx, -margy, pos(3)+2*margx, pos(4)+2*margy];
num_getframes = 5000;
frames_new = cell(num_getframes,1);
for n_fr = 1%:num_getframes%num_frames2
    %pause(1/60);
    im1.CData = frames_cut_df(:,:,n_fr);
    im_yup.YData = yup_loc_f(n_fr);
    im_ydn.YData = ydn_loc_f(n_fr);
    drawnow();
    frames_new{n_fr} = getframe(f1,rect1);
end

% figure; imagesc(frames_new{1}.cdata); axis off; 
% 
% writerObj = VideoWriter([fpath, 'pupil_track_out.avi']); 
% writerObj.FrameRate = 60;
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for n_fr=1:length(frames_new)
%     % convert the image to a frame
%     frame = frames_new{n_fr}.cdata;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);
%%
clear frames_cut_df
%%

video_input2 = VideoReader([fpath, 'A2_ammn1_body0001 21-04-11 13-13-05.avi']);

%%

num_frames = video_input2.NumFrames;
%num_frames = 5000;

%%

frame1 = rgb2gray(video_input2.readFrame);

fr2 = video_input2.FrameRate;

[d1, d2] = size(frame1);

frames = zeros([d1, d2, num_frames], 'uint8');

frames(:,:,1) = frame1;

for n_fr = 2:num_frames
    frames(:,:,n_fr) = rgb2gray(video_input2.readFrame);
end

%%

figure; imagesc(mean(frames(:,:,1:500),3)); axis equal tight;
title('Mean video; laser region (2 clicks)')
[x_laser, y_laser] = ginput(2);
x_laser = round(x_laser);
y_laser = round(y_laser);
rectangle('Position',[min(x_laser) min(y_laser) diff(x_laser)  diff(y_laser)],'EdgeColor', 'r')
title('Mean video; whisker region (2 clicks)')
[x_wh, y_wh] = ginput(2);
x_wh = round(x_wh);
y_wh = round(y_wh);
rectangle('Position',[min(x_wh) min(y_wh) diff(x_wh)  diff(y_wh)],'EdgeColor', 'g')
title('Mean video; paw region (2 clicks)')
[x_paw, y_paw] = ginput(2);
x_paw = round(x_paw);
y_paw = round(y_paw);
rectangle('Position',[min(x_paw) min(y_paw) diff(x_paw)  diff(y_paw)],'EdgeColor', 'y')
title('Mean video full')

%%
mean_laser_int = squeeze(mean(mean(frames(min(y_laser):min(y_laser)+diff(y_laser),min(x_laser):min(x_laser)+diff(x_laser),:))));
whisker_frames = frames(min(y_wh):min(y_wh)+diff(y_wh),min(x_wh):min(x_wh)+diff(x_wh),:);
paw_frames = frames(min(y_paw):min(y_paw)+diff(y_paw),min(x_paw):min(x_paw)+diff(x_paw),:);


%%
figure; plot(mean_laser_int)
title('Mean laser activity; Select laser on thresh (1 click)')
[~, thresh1] = ginput(1);
title('Mean laser activity;')

temp_vals = find(mean_laser_int>thresh1);
onset_fr = temp_vals(1);

whisker_frames2 = whisker_frames(:,:,onset_fr:end);
paw_frames2 = paw_frames(:,:,onset_fr:end);

clear whisker_frames paw_frames;

%% do pca of whisker and paw frames to get some motion energy

[d1, d2, num_frames_w] = size(whisker_frames2);

whisker_frames_2d = single(reshape(whisker_frames2, [d1*d2, num_frames_w]));

%%

[wh_coeff,wh_score,wh_latent,wh_tsquared,wh_explained,wh_mu] = pca(whisker_frames_2d);

%%
lr_pcs = 1:3;
figure; plot(wh_coeff(:,lr_pcs))

%%
mov_low_rank = reshape(wh_score(:,lr_pcs)*wh_coeff(1:1000,lr_pcs)',d1,d2,[]);
f_save_tif_stack2_YS(mov_low_rank, [fpath, sprintf('whisker_comp %s', strjoin(string(lr_pcs)))])

%%

x_w = (1:num_frames_w)/fr2;
xq = proc_data.data.frame_data.frame_times_mpl{1,1}/1000;

w_coeff2_itp = interp1(x_w,wh_coeff(:,2),xq);
w_coeff3_itp = interp1(x_w,wh_coeff(:,3),xq);

%%

base_onset_win = [5 35];
%tt = 270;
tt = [170, 270, 3, 6];

col1 = {[1 .3 .3], [1 .3 .3], [.7 .7 .7], [.7 .7 .7]};
col2 = {'r', 'r', 'k', 'k'};
%%
stim_frame_index = proc_data.data.stim_frame_index{1};
trial_types = proc_data.data.trial_types;
MMN_orientations = proc_data.data.MMN_orientations;
w_coeff3_itp_sm = f_smooth_gauss2(w_coeff3_itp', 3, 0)';
figure; plot(xq,w_coeff3_itp); hold on;
plot(xq,w_coeff3_itp_sm)
xlabel('sec')
t_p = (-4:base_onset_win(2))/10;

figure; hold on; axis tight;
wh_mean_tr = cell(4,1);
wh_trial_means = cell(4,1);
for n_tt = 1:4
    tt_temp = tt(n_tt);
    stim_frame_index1 = stim_frame_index(trial_types == tt_temp);
    num_trials = numel(stim_frame_index1);
    
    
    
    w_coeff_sort = squeeze(f_get_stim_trig_resp(w_coeff3_itp_sm', stim_frame_index1, base_onset_win));
    w_coeff_sort_n = w_coeff_sort - mean(w_coeff_sort(1:5,:));
    
    
    
%     rem_thr = -80;
%     rem_trial = false(num_trials,1);
%     for n_tr = 1:num_trials
%         if sum(pupil_sort_n(:,n_tr)<rem_thr)
%             rem_trial(n_tr) =1;
%         end
%     end
    
    %pupil_sort_n2 = pupil_sort_n(:,~rem_trial);
    %pupil_sort2 = pupil_sort(:,~rem_trial);
    
    plot(t_p, w_coeff_sort, '--', 'color', col1{n_tt})
    wh_mean_tr{n_tt} = mean(w_coeff_sort,2);
    
    wh_trial_means{n_tt} = mean(w_coeff_sort);
end
for n_tt = 1:4
    plot(t_p, wh_mean_tr{n_tt}, 'Color', col2{n_tt}, 'Linewidth', 2)
end
xlabel('Time (sec)'); ylabel('smoothed whisking mag')

figure; hold on;
for n_tt = 1:4
    mean2 = mean(pup_mean_tr{n_tt});
    plot(n_tt, pup_mean_tr{n_tt}, 'o', 'Color', col2{n_tt})
    
end


figure; plot(w_coeff3_itp,pupil_diam_ds, '.')
xlabel('whisking frames'); ylabel('pupil frames')