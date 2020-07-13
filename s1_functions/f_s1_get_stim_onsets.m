function stim_times_trace = f_s1_get_stim_onsets(stim_volt_trace, thresh)

num_samp = numel(stim_volt_trace);
stim_times_trace = zeros(size(stim_volt_trace));
for n_samp = 2:num_samp
    if (stim_volt_trace(n_samp) > thresh) && (stim_volt_trace(n_samp-1) < thresh)
        stim_times_trace(n_samp) = 1;
    end
end

% figure;  hold on;
% plot(stim_volt_trace)
% plot(stim_times_trace)

end