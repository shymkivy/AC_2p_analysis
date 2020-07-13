function shift = f_s1_align_traces_xcor(ca_traces, frame_times, volt_data)

v_times = 1:size(volt_data,1);

% figure; hold on;
% plot(v_times, volt_data);
% plot(frame_times,ca_traces);

% upsample ca data
x = frame_times;
v = ca_traces;
xq = v_times;
vq = interp1(x, v, xq);
vq(isnan(vq)) = 0;
% figure; hold on; 
% plot(xq,vq); 
% plot(frame_times,ca_traces);

[xcor1, lags] = xcorr(vq, volt_data);
[~, ind_max] = max(xcor1);
shift = lags(ind_max);

% figure;
% plot(xq+shift,volt_data);
% hold on;
% plot(frame_times, ca_traces);

end