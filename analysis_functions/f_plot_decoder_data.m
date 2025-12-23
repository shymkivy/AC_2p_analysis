function f_plot_decoder_data(dec_data, plot_t)
% num_dsets, num_frames, num_trials

dec_acc_frames = cat(1,dec_data.accuracy);
dec_acc_frames_shuff = cat(1,dec_data.accuracy_shuff);

figure; hold on; axis tight
plot(plot_t, dec_acc_frames_shuff', color=[0 0 0 0.2])
plot(plot_t, mean(dec_acc_frames_shuff,1), color=[0 0 0], LineWidth=2)

plot(plot_t, dec_acc_frames', color=[0 0.4470 0.7410 0.2])
plot(plot_t, mean(dec_acc_frames,1), color=[0 0.4470 0.7410], LineWidth=2)
title("Binwise decoder");
ylabel("Performance");
xlabel("Time");
ylim([0, 1]);
end