clear;
close all;

fpath = 'C:\Users\ys2605\Desktop\stuff\AC_data\1_30_21\stim_fov2_cell4-011';
fname = 'stim_fov2_cell4-011_Cycle00001_VoltageRecording_001.csv';
fname_roi = 'cell_roi_trace.csv';

framePeriod = 0.01744251 * 1000;

stim_chan = 5;

%%
addpath([pwd '\general_functions']);
%%
data = csvread([fpath '\' fname], 1, 0);
data_stim = csvread([fpath '\' fname_roi], 1, 0);

roi_data_norm = (data_stim(:,2) - min(data_stim(:,2)))/max(data_stim(:,2) - min(data_stim(:,2)));
roi_data_norm_sm = smooth(roi_data_norm, 10);

stim_trace_norm = data(:,stim_chan+1)/max(data(:,stim_chan+1));

figure; hold on; axis tight
plot(data(:,1), stim_trace_norm);
plot(data_stim(:,1)*framePeriod, roi_data_norm);
plot(data_stim(:,1)*framePeriod, roi_data_norm_sm, 'k', 'LineWidth', 1);


stim_trace_norm_dsmp = interp1(data(:,1),stim_trace_norm, data_stim(:,1)*framePeriod);
onset_frames = f_get_pulse_times(stim_trace_norm_dsmp,0.5, 999);

stim_t = ((1:150)-30)*framePeriod/1000;
resp_sort = squeeze(f_get_stim_trig_resp(roi_data_norm_sm', onset_frames, [30 120]));
stim_sort = squeeze(f_get_stim_trig_resp(stim_trace_norm_dsmp', onset_frames, [30 120]));

resp_sort_bs = resp_sort - mean(resp_sort(1:30,:));
peak_val = max(resp_sort_bs(:));

figure;
imagesc([resp_sort_bs, stim_sort(:,1)*peak_val]');
title('stim trig response');



figure; hold on; axis tight;
plot(stim_t, resp_sort_bs, 'color', [.5 .5 .5])
plot(stim_t, mean(resp_sort_bs,2), 'm', 'LineWidth', 2);
plot(stim_t, stim_sort*peak_val, 'r');
title('stim trig response');
legend('trials', 'trial ave', 'stim');