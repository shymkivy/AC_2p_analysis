function f_plot_ctx_traces2(data, ops, area1, dataset_num)

trial_ave_z = data.proc{dataset_num}.trial_ave_z;
context_resp_cells = data.proc{dataset_num}.context_resp_cells;
MMN_freq = data.stim_params{dataset_num}.MMN_freq;

MMN_freq = MMN_freq(1,:);

plot_trace = trial_ave_z(:,ops.win.analysis_window,:); %ops.time_stim_window>=0
plot_time = ops.time_stim_window(ops.win.analysis_window);
plot_contx_traces(plot_time, plot_trace,[], MMN_freq, [0 max(max(mean(plot_trace,1)+0.5))], ops)
suptitle(sprintf('%s, dataset %d', area1, dataset_num));

plot_trace = trial_ave_z(:,ops.win.analysis_window,:); %ops.time_stim_window>=0
plot_time = ops.time_stim_window(ops.win.analysis_window);
plot_contx_traces(plot_time, plot_trace, context_resp_cells, MMN_freq, [0 max(max(mean(plot_trace,1)+0.5))], ops)
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

function plot_contx_traces(plot_time, plot_trace, ctx_cells, MMN_freq, ylims1, ops)
figure;

% first line
for ii = 1:10
    subplot(3,10,ii)
    cute_plot(plot_time, mean(plot_trace(:,:,ii),1), ops);
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
if isempty(ctx_cells)
    cute_plot(plot_time, mean(plot_trace(:,:,MMN_freq(2)),1), ops);
else
    cute_plot(plot_time, mean(plot_trace(ctx_cells(:,1),:,MMN_freq(2)),1), ops);
end
ylim(ylims1);
title(sprintf('cont'));
ylabel(sprintf('freq %d', MMN_freq(2)))
ax = gca;
ax.XTick = temp_tick;
for ii = 11:18
    subplot(3,10,ii+1);
    if isempty(ctx_cells)
        cute_plot(plot_time, mean(plot_trace(:,:,ii),1), ops);
    else
        cute_plot(plot_time, mean(plot_trace(ctx_cells(:,1),:,ii),1), ops);
    end
    ylim(ylims1);
    title(sprintf('red %d', ii-10));
    ax = gca;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
end
subplot(3,10,20);
if isempty(ctx_cells)
    cute_plot(plot_time, mean(plot_trace(:,:,19),1), ops);
else
    cute_plot(plot_time, mean(plot_trace(ctx_cells(:,1),:,19),1), ops);
end
ylim(ylims1);
title(sprintf('dev'));
ax = gca;
ax.YTickLabel = [];
ax.XTickLabel = [];
% line 3
subplot(3,10,21);
if isempty(ctx_cells)
    cute_plot(plot_time, mean(plot_trace(:,:,MMN_freq(1)),1), ops);
else
    cute_plot(plot_time, mean(plot_trace(ctx_cells(:,2),:,MMN_freq(1)),1), ops);
end
ylim(ylims1);
title(sprintf('cont'));
ylabel(sprintf('freq %d', MMN_freq(1)))
ax = gca;
ax.XTick = temp_tick;
for ii = 20:27
    subplot(3,10,ii+2);
    if isempty(ctx_cells)
        cute_plot(plot_time, mean(plot_trace(:,:,ii),1), ops);
    else
        cute_plot(plot_time, mean(plot_trace(ctx_cells(:,2),:,ii),1), ops);
    end
    ylim(ylims1);
    title(sprintf('red %d', ii-19));
    ax = gca;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
end
subplot(3,10,30);
if isempty(ctx_cells)
    cute_plot(plot_time, mean(plot_trace(:,:,28),1), ops);
else
    cute_plot(plot_time, mean(plot_trace(ctx_cells(:,2),:,28),1), ops);
end
ylim(ylims1);
title(sprintf('dev'));
ax = gca;
ax.YTickLabel = [];
ax.XTickLabel = [];

end

function cute_plot(plot_time, plot_trace, ops)
    plot(plot_time, plot_trace, 'Color', [1 0 1], 'LineWidth', 2);
end
