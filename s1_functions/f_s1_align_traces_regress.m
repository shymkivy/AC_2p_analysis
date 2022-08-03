function [shift, scaling_factor] = f_s1_align_traces_regress(ca_traces, frame_times, volt_data, thresh)

%v_times = 1:size(volt_data,1);
% figure; hold on;
% plot(v_times, volt_data);
% plot(frame_times,ca_traces);
% 
if ~exist('thresh', 'var')
    thresh = [.3 .8];
end
%%
on_times_low = [];
on_times_high = [];
for n_fr = 2:numel(ca_traces)
    if and(ca_traces(n_fr) > thresh(1), ca_traces(n_fr-1) < thresh(1))
        on_times_low = [on_times_low; n_fr-1];
    end
    if and(ca_traces(n_fr) > thresh(2), ca_traces(n_fr-1) < thresh(2))
        on_times_high = [on_times_high; n_fr];
    end
end

while numel(on_times_low) ~= numel(on_times_high)
    figure; plot(ca_traces); axis tight;
    title(sprintf('threshhold [%.1f %.1f] failed, select manual (2 clicks)', thresh(1), thresh(2)));
    [~,thresh] = ginput(2);
    close;
    
    on_times_low = [];
    on_times_high = [];
    for n_fr = 2:numel(ca_traces)
        if and(ca_traces(n_fr) > thresh(1), ca_traces(n_fr-1) < thresh(1))
            on_times_low = [on_times_low; n_fr-1];
        end
        if and(ca_traces(n_fr) > thresh(2), ca_traces(n_fr-1) < thresh(2))
            on_times_high = [on_times_high; n_fr];
        end
    end
end

pulse_times_on = (frame_times(on_times_high) + frame_times(on_times_low))/2 - mean(diff(frame_times))/2;

%%
off_times_low = [];
off_times_high = [];
for n_fr = 2:numel(ca_traces)
    if and(ca_traces(n_fr) < thresh(2), ca_traces(n_fr-1) > thresh(2))
        off_times_low = [off_times_low; n_fr-1];
    end
    if and(ca_traces(n_fr) < thresh(1), ca_traces(n_fr-1) > thresh(1))
        off_times_high = [off_times_high; n_fr];
    end
end
pulse_times_off = (frame_times(off_times_high) + frame_times(off_times_low))/2 - mean(diff(frame_times))/2;

%%
on_times_low = [];
on_times_high = [];
for n_t = 2:numel(volt_data)
    if and(volt_data(n_t) > thresh(1), volt_data(n_t-1) < thresh(1))
        on_times_low = [on_times_low; n_t-1];
    end
    if and(volt_data(n_t) > thresh(2), volt_data(n_t-1) < thresh(2))
        on_times_high = [on_times_high; n_t];
    end
end
pulse_times_on_volt = (on_times_high + on_times_low)/2 - .5;

%%
off_times_low = [];
off_times_high = [];
for n_t = 2:numel(volt_data)
    if and(volt_data(n_t) < thresh(2), volt_data(n_t-1) > thresh(2))
        off_times_low = [off_times_low; n_t-1];
    end
    if and(volt_data(n_t) < thresh(1), volt_data(n_t-1) > thresh(1))
        off_times_high = [off_times_high; n_t];
    end
end
pulse_times_off_volt = (off_times_high + off_times_low)/2 - .5;

%%
%x_on = [ones(size(pulse_times_on_volt)), pulse_times_on_volt]\pulse_times_on;
%x_off = [ones(size(pulse_times_off_volt)), pulse_times_off_volt]\pulse_times_off;

pulse_comb = [pulse_times_on; pulse_times_off];
pulse_comb_volt = [pulse_times_on_volt; pulse_times_off_volt];
x_comb = [ones(size(pulse_comb_volt)), pulse_comb_volt]\pulse_comb;

%%
scaling_factor = x_comb(2);
shift = x_comb(1);

%%

% figure; hold on;
% plot(frame_times, ca_traces)
% plot(pulse_times_on, zeros(size(pulse_times_on)), 'o')
% plot(pulse_times_off, ones(size(pulse_times_on)), 'o')
% plot((1:numel(volt_data))*x_comb(2)+x_comb(1), volt_data) % 

end
