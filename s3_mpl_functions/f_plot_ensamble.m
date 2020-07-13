function f_plot_ensamble(ens1, trial_data_sort_act, trial_types_dred, n_comp1)

tt1 = unique(trial_types_dred);
figure;
sp = cell(numel(tt1),1);
ymin = 0; ymax = 0.01;
for n_trial = 1:numel(tt1)
    sp{n_trial} = subplot(2,ceil(numel(tt1)/2),n_trial);
    plot_trace = mean(mean(trial_data_sort_act(ens1,:,trial_types_dred==tt1(n_trial)),3),1)';
    plot(plot_trace);
    ymin = min([ymin; plot_trace(:)]);
    ymax = max([ymax; plot_trace(:)]);
    title(num2str(tt1(n_trial)))
end
for n_trial = 1:numel(tt1)
    ylim(sp{n_trial}, [ymin ymax]);
end
suptitle(sprintf('comp %d, %d cells', n_comp1, sum(ens1)))


end