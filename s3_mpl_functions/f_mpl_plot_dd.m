function [plot_hand,ylim1]  = f_mpl_plot_dd(trial_ave, resp_cells, trial_t, ops, plot_hand)
% input has to be 'c rf d c r df' or  c1 203 170 c2 103 270

if ~exist('plot_hand', 'var') || ~iscell(plot_hand)
    if ops.plot_combined
        m_pl = 3;
    else
        m_pl = 2;
    end
    plot_hand = cell(m_pl,1);
    figure;
    for n_flip = 1:m_pl
        plot_hand{n_flip} = subplot(1,m_pl,n_flip); hold on;
    end
end


%%
resp_cells1 = reshape(resp_cells, [], 3,2);

traces1 = cell(m_pl,1);
traces1{1} = trial_ave(logical(sum(resp_cells1(:,:,1),2)),:,[1 2 3]);
traces1{2} = trial_ave(logical(sum(resp_cells1(:,:,2),2)),:,[4 5 6]);
if ops.plot_combined
    traces1{3} = cat(1, traces1{1}, traces1{2});
end

ymin = 0;
ymax = 0.02;
title1 = {'c rf d', 'c r df', 'combined'};
for n_flip = 1:m_pl
    traces2 = traces1{n_flip};
    num_cells = size(traces2,1);
    subplot(1,m_pl,n_flip); hold on;
    for n_ctx = 1:3
        temp_trace = traces2(:,:,n_ctx);
        temp_trace_mean = mean(temp_trace);
        temp_trace_sem = std(temp_trace)/sqrt(num_cells-1);
        shadedErrorBar_YS(trial_t, temp_trace_mean, temp_trace_sem, ops.context_colors{n_ctx});
        ymin = min([ymin (temp_trace_mean-temp_trace_sem)]);
        ymax = max([ymax (temp_trace_mean+temp_trace_sem)]);
    end
    axis tight;
    title([title1{n_flip} '; ' num2str(num_cells) ' cells']);
end

ylim1 = [ymin ymax];
for n_flip = 1:numel(plot_hand)
    ylim(plot_hand{n_flip}, ylim1);
end

end