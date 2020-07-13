function f_plot_trials(trial_ave, trials_to_plot, trial_t, plot_hand)

num_trials = numel(trials_to_plot);
num_cells = size(trial_ave,1);

if ~exist('plot_hand', 'var') || ~iscell(plot_hand)
    plot_hand = cell(num_trials,1);
    figure;
    [m,n] = f_give_subplotdims(num_trials);
    for n_tr = 1:num_trials
        plot_hand{n_tr} = subplot(m,n,n_tr); hold on;
    end
end




%title1 = {'c rf d', 'c r df', 'combined'};
ymin = 0;
ymax = 0.02;
for n_tr = 1:num_trials
    subplot(plot_hand{n_tr});
    temp_ave = squeeze(trial_ave(:,:,n_tr));
    temp_sem = std(temp_ave)/sqrt(num_cells-1);
    shadedErrorBar(trial_t, mean(temp_ave), temp_sem); %, 'lineprops', ops.context_colors{n_ctx});
    ymin = min([ymin (mean(temp_ave) - temp_sem)]);
    ymax = max([ymax (mean(temp_ave) + temp_sem)]);
    axis tight;
    title(['Freq' num2str(trials_to_plot(n_tr))])
end
for n_tr = 1:num_trials
    ylim(plot_hand{n_tr}, [ymin ymax]);
end

end