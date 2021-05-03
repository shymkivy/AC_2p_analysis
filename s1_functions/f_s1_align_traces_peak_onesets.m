function [shift, scaling_factor] = f_s1_align_traces_peak_onesets(ca_traces, frame_times, volt_data)

%v_times = 1:size(volt_data,1);
% figure; hold on;
% plot(v_times, volt_data);
% plot(frame_times,ca_traces);
% 

thresh = 0.5;

[pulse_on_ca_idx, pulse_off_ca_idx] = f_get_pulse_times(ca_traces,thresh);
pulse_on_ca = frame_times(pulse_on_ca_idx-1);
pulse_off_ca = frame_times(pulse_off_ca_idx);

[pulse_on_volt, pulse_off_volt] = f_get_pulse_times(volt_data,thresh);
pulse_on_volt = pulse_on_volt - 1;

%%
pulse_on_ca2 = pulse_on_ca - pulse_on_ca(1);
pulse_off_ca2 = pulse_off_ca - pulse_off_ca(1);
pulse_ca_mean = mean([pulse_on_ca2, pulse_off_ca2],2);

pulse_on_volt2 = pulse_on_volt - pulse_on_volt(1);
pulse_off_volt2 = pulse_off_volt - pulse_off_volt(1);
pulse_volt_mean = mean([pulse_on_volt2, pulse_off_volt2],2);

scaling_factor = mean(pulse_ca_mean(2:end)./pulse_volt_mean(2:end));
%%
shift = mean(pulse_on_ca - pulse_on_volt*scaling_factor);

end
