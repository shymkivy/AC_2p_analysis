function f_mpl_plot_ctx(trial_ave, MMN_freq,plot_time,ops)

plot_trace = squeeze(mean(trial_ave,1)); 

if size(MMN_freq,1) > 1
    [num_cells, num_bins, ~] = size(trial_ave);
    cont_plot = zeros(num_cells,num_bins,2);
    for n_cell = 1:num_cells
        cont_plot(n_cell,:,1) = trial_ave(n_cell,:,MMN_freq(n_cell,2));
        cont_plot(n_cell,:,2) = trial_ave(n_cell,:,MMN_freq(n_cell,1));
    end
    mmn_plot1 = squeeze(mean(cont_plot(:,:,1)));
    mmn_plot2 = squeeze(mean(cont_plot(:,:,2)));
else
    mmn_plot1 = plot_trace(:,MMN_freq(2));
    mmn_plot2 = plot_trace(:,MMN_freq(1));
end

ylims1 = [min([plot_trace(:); mmn_plot1(:); mmn_plot2(:)]) max([plot_trace(:); mmn_plot1(:); mmn_plot2(:)])];

figure;
% first line
for ii = 1:10
    subplot(3,10,ii)
    cute_plot(plot_time, plot_trace(:,ii));
    ylim(ylims1);
    title(sprintf('Freq %d', ii));
    ax = gca;
    if ii == 1
        temp_tick = ax.XTick;
    elseif ii > 1
%         ax.YTickLabel = [];
%         ax.XTickLabel = [];
    end
end
% second line
for ii = 11:18
    subplot(3,10,ii);
    cute_plot(plot_time, plot_trace(:,ii));
    ylim(ylims1);
    title(sprintf('RedF %d', ii-10));
%    ax = gca;
%     ax.YTickLabel = [];
%     ax.XTickLabel = [];
end
subplot(3,10,19);
cute_plot(plot_time, plot_trace(:,20));
ylim(ylims1);
title(sprintf('Dev'));
%ax = gca;
% ax.YTickLabel = [];
% ax.XTickLabel = [];



subplot(3,10,20); hold on;
plot(plot_time, plot_trace(:,19), 'Color', 'b', 'LineWidth', 2);
plot(plot_time, plot_trace(:,20), 'Color', 'r', 'LineWidth', 2);
plot(plot_time, mmn_plot1, 'Color', 'k', 'LineWidth', 2);
axis tight;
ylim(ylims1);
title(sprintf('c fr d'));
%ax = gca;
%ax.XTick = temp_tick;

% line 3

for ii = 21:28
    subplot(3,10,ii);
    cute_plot(plot_time, plot_trace(:,ii));
    ylim(ylims1);
    title(sprintf('Red %d', ii-20));
%    ax = gca;
%     ax.YTickLabel = [];
%     ax.XTickLabel = [];
end

subplot(3,10,29);
cute_plot(plot_time, plot_trace(:,30));
ylim(ylims1);
title(sprintf('DevF'));
%ax = gca;
% ax.YTickLabel = [];
% ax.XTickLabel = [];

subplot(3,10,30); hold on;
plot(plot_time, plot_trace(:,29), 'Color', 'b', 'LineWidth', 2);
plot(plot_time, plot_trace(:,30), 'Color', 'r', 'LineWidth', 2);
plot(plot_time, mmn_plot2, 'Color', 'k', 'LineWidth', 2);
axis tight;
ylim(ylims1);
title(sprintf('c r fd'));
%ax = gca;
% ax.YTickLabel = [];
% ax.XTickLabel = [];



end

function cute_plot(plot_time, plot_trace)
    plot(plot_time, plot_trace, 'Color', [1 0 1], 'LineWidth', 2);
    axis tight;
end
