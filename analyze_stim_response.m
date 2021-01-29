
fpath = 'C:\Users\ys2605\Desktop\stuff\AC_data\1_3_21_stim\rest_fov3_cell8-010';
fname = 'rest_fov3_ce8-010_Cycle00001_VoltageRecording_001.csv';
fname_roi = 'cell_roi_trace.csv';

framePeriod = 0.01744251 * 1000;

stim_chan = 1;

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
resp_sort = f_get_stim_trig_resp(roi_data_norm_sm', onset_frames, [30 120]);
stim_sort = f_get_stim_trig_resp(stim_trace_norm_dsmp', onset_frames, [30 120]);
figure;
imagesc([squeeze(resp_sort)'; squeeze(stim_sort(1,:,1))]);
title('stim trig response');


figure; hold on; axis tight;
plot(stim_t, squeeze(resp_sort), 'k');
plot(stim_t, mean(squeeze(resp_sort),2), 'm', 'LineWidth', 2);
plot(stim_t, squeeze(stim_sort), 'r');
title('stim trig response');
legend('trials', 'trial ave', 'stim');