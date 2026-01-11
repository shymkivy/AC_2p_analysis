function [trial_window_t, num_baseline_resp_frames] = f_dv_compute_window_t(trial_window, vol_period)
vol_period2 = vol_period/1000;

frame_start = ceil(trial_window(1)/vol_period2);
frame_end = floor(trial_window(2)/vol_period2);
trial_window_t = (frame_start:frame_end)*vol_period2;
num_baseline_resp_frames = [frame_start frame_end];     

end