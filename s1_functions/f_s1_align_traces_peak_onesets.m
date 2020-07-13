function shift = f_s1_align_traces_peak_onesets(ca_traces, frame_times, volt_data)

%v_times = 1:size(volt_data,1);
% figure; hold on;
% plot(v_times, volt_data);
% plot(frame_times,ca_traces);
% 

thresh = 0.5;

[pulse_on_ca_idx, pulse_off_ca_idx] = f_get_pulse_times(ca_traces,thresh);
pulse_on_ca = frame_times(pulse_on_ca_idx);
pulse_off_ca = frame_times(pulse_off_ca_idx);

[pulse_on_volt, pulse_off_volt] = if_get_pulse_times(volt_data,thresh);


onset_shift = mean(pulse_on_ca - pulse_on_volt);
offset_shift = mean(pulse_off_ca - pulse_off_volt);

shift = mean([onset_shift,offset_shift]);

end
