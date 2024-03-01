function f_plot_raster_mean(raster, plot_mean, xlabels, trial_ind, colors1, color_map, invert_cmap)

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

if ~exist('colors1', 'var')
    colors1 = [];
end

if ~exist('xlabels', 'var')
    xlabels = [];
end

if ~exist('color_map', 'var')
    color_map = 'parula';
end

if ~exist('invert_cmap', 'var')
    invert_cmap = 0;
end

if invert_cmap
    raster2 = -raster;
else
    raster2 = raster;
end

raster2n = raster2-min(raster2(:));
raster2n = raster2n./max(raster2n(:));

figure;
if ~plot_mean
    imagesc(xlabels, [], raster2n);
    axis tight;
    ylabel('Cells');
else
    s1 = subplot(4,1,1:3);
    
    imagesc(xlabels, [], raster2n);
    axis tight;
    ylabel('Cells');
    s1.XAxis.TickValues = [];
    colormap(color_map)
    if plot_tr_ind
        f_plot_trial_indicator3(raster2n, trial_ind, 1, colors1, xlabels)
    end
    s2 = subplot(4,1,4);
    plot(xlabels, mean(raster), 'k');
    ylabel('Population average');
    xlabel('Time (s)')
    %s2.XAxis.TickValues = [];
    linkaxes([s1,s2],'x')
    subplot(s1); axis tight;
end


end