function f_plot_ctx_traces(data, ops, area1, dataset_num)

cells_ind = data.(area1).num_dataset ==dataset_num;
trial_ave = data.(area1).trial_ave_z(cells_ind,:,:);
ctx_tuned_cells = logical(sum(data.(area1).context_resp_cells(cells_ind,:),2));


MMN_freq = data.(area1).MMN_freq(cells_ind,:);
MMN_freq = MMN_freq(1,:);


plot_trace = mean(trial_ave(:,ops.win.analysis_window,:),1); %ops.time_stim_window>=0
plot_time = ops.time_stim_window(ops.win.analysis_window);
plot_contx_traces(plot_time, plot_trace,MMN_freq, [0 max(plot_trace(:)+0.5)], ops)
suptitle(sprintf('%s, dataset %d', area1, dataset_num));

plot_trace = mean(trial_ave(ctx_tuned_cells,ops.win.analysis_window,:),1); %ops.time_stim_window>=0
plot_time = ops.time_stim_window(ops.win.analysis_window);
plot_contx_traces(plot_time, plot_trace,MMN_freq, [0 max(plot_trace(:)+0.5)], ops)
suptitle(sprintf('%s, ctx dataset %d', area1, dataset_num));

% cells_loc = data.A1.cell_loc(cells_ind);
% cell_im = zeros(256,256);
% for n_cell = 1:numel(cells_loc)
%     temp_coord = cells_loc{n_cell};
%     for jj = 1:size(temp_coord,1)
%         cell_im(temp_coord(jj,1),temp_coord(jj,2)) = 1;
%     end
% end
% 
% figure;
% imagesc(cell_im)

end

function plot_contx_traces(plot_time, plot_trace, MMN_freq, ylims1, ops)
figure;

% first line
for ii = 1:10
    subplot(3,10,ii)
    cute_plot(plot_time, plot_trace(:,:,ii), ops);
    ylim(ylims1);
    title(sprintf('freq %d', ii));
    ax = gca;
    if ii == 1
        temp_tick = ax.XTick;
    elseif ii > 1
        ax.YTickLabel = [];
        ax.XTickLabel = [];
    end
end
% second line
subplot(3,10,11);
cute_plot(plot_time, plot_trace(:,:,MMN_freq(2)), ops);
ylim(ylims1);
title(sprintf('cont'));
ylabel(sprintf('freq %d', MMN_freq(2)))
ax = gca;
ax.XTick = temp_tick;
for ii = 11:18
    subplot(3,10,ii+1);
    cute_plot(plot_time, plot_trace(:,:,ii), ops);
    ylim(ylims1);
    title(sprintf('red %d', ii-10));
    ax = gca;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
end
subplot(3,10,20);
cute_plot(plot_time, plot_trace(:,:,19), ops);
ylim(ylims1);
title(sprintf('dev'));
ax = gca;
ax.YTickLabel = [];
ax.XTickLabel = [];
% line 3
subplot(3,10,21);
cute_plot(plot_time, plot_trace(:,:,MMN_freq(1)), ops);
ylim(ylims1);
title(sprintf('cont'));
ylabel(sprintf('freq %d', MMN_freq(1)))
ax = gca;
ax.XTick = temp_tick;
for ii = 20:27
    subplot(3,10,ii+2);
    cute_plot(plot_time, plot_trace(:,:,ii), ops);
    ylim(ylims1);
    title(sprintf('red %d', ii-19));
    ax = gca;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
end
subplot(3,10,30);
cute_plot(plot_time, plot_trace(:,:,28), ops);
ylim(ylims1);
title(sprintf('dev'));
ax = gca;
ax.YTickLabel = [];
ax.XTickLabel = [];

end

function cute_plot(plot_time, plot_trace, ops)
    plot(plot_time, plot_trace, 'Color', [1 0 1], 'LineWidth', 2);
end
