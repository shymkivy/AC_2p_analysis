function f_plot_raster_mean(raster, cell_seq, trial_seq, ops, plot_mean, trial_ind)


[num_cells, num_trials] = size(raster);

if isempty(cell_seq)
    cell_seq = 1:num_cells;
end
if isempty(trial_seq)
    trial_seq = 1:num_trials;
end
if ~exist('trial_ind', 'var') || isempty(trial_ind)
    plot_tr_ind = 0;
else
    plot_tr_ind = 1;
end
if ~exist('plot_mean', 'var') || isempty(plot_mean)
    plot_mean = 0;
else
    plot_mean = 1;
end

raster2 = raster(cell_seq,trial_seq);

figure;
if ~plot_mean
    imagesc(raster2);
else
    s1 = subplot(4,1,1:3);
    imagesc(raster2); axis tight;
    if plot_tr_ind
        f_plot_trial_indicator3(raster2, trial_ind, 1, ops)
    end
    s2 = subplot(4,1,4);
    plot(mean(raster2), 'k');
    s2.XAxis.TickValues = []
    linkaxes([s1,s2],'x')
    subplot(s1); axis tight;
end


end