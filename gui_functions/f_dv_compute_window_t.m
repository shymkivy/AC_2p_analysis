function [trial_window_t, num_baseline_resp_frames] = f_dv_compute_window_t(trial_window, vol_period)
vol_period2 = vol_period/1000;

trial_window_t = (ceil(trial_window(1)/vol_period2):floor(trial_window(2)/vol_period2))*vol_period2;
num_baseline_resp_frames = [sum(trial_window_t<=0) sum(trial_window_t>0)];     

end