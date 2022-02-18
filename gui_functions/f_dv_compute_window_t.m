function [trial_window_t, num_baseline_resp_frames] = f_dv_compute_window_t(app, trial_window)

frame_period = 1/app.FramerateEditField.Value;

trial_window_t = (ceil(trial_window(1)/frame_period):floor(trial_window(2)/frame_period))*frame_period;
num_baseline_resp_frames = [sum(trial_window_t<=0) sum(trial_window_t>0)];     


end