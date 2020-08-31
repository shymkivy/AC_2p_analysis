function f_mpl_plot_ctx2(trial_ave, resp_cells, ctx_mmn,plot_time,ops)

if ops.plot_combined
    m_pl = 3;
else
    m_pl = 2;
end

[num_cells, num_bins, ~] = size(trial_ave);
col_list = [ops.context_colors{1}; repmat(ops.context_colors{2}, 8, 1); ops.context_colors{3}];
resp_cells1 = reshape(resp_cells, [], 3,2);
leg_list = cat(1,sprintf('cont %d', ctx_mmn(1)), ops.context_types_labels(11:18), ops.context_types_labels(20),...
            sprintf('cont %d', ctx_mmn(4)), ops.context_types_labels(21:28), ops.context_types_labels(30));

traces1 = cell(m_pl,1);
if size(ctx_mmn,1) > 1
    traces1{1} = zeros(num_cells,num_bins,10);
    traces1{2} = zeros(num_cells,num_bins,10);
    for n_cell = 1:num_cells
        traces1{1}(n_cell,:,:) = trial_ave(n_cell,:,[ctx_mmn(n_cell,1), 11:18, ctx_mmn(n_cell,3)]);
        traces1{2}(n_cell,:,:) = trial_ave(n_cell,:,[ctx_mmn(n_cell,4), 21:28, ctx_mmn(n_cell,6)]);
    end
else
    traces1{1} = trial_ave(:,:,[ctx_mmn(1), 11:18, ctx_mmn(3)]);
    traces1{2} = trial_ave(:,:,[ctx_mmn(4), 21:28, ctx_mmn(6)]);
end
traces1{1} = traces1{1}(logical(sum(resp_cells1(:,:,1),2)),:,:);
traces1{2} = traces1{2}(logical(sum(resp_cells1(:,:,2),2)),:,:);
if ops.plot_combined
    traces1{3} = cat(1, traces1{1}, traces1{2});
end

plot_hand = cell(m_pl*10,1);
ylims1 = [0 0.2];

figure;
for n_flip = 1:m_pl
    traces2 = traces1{n_flip};
    num_cells2 = size(traces2,1);
    for n_tr = 1:10
        n_pl = 10*(n_flip-1)+n_tr;
        plot_hand{n_pl} = subplot(m_pl,10,n_pl);
        temp_trace = traces2(:,:,n_tr);
        temp_trace_mean = mean(temp_trace);
        temp_trace_sem = std(temp_trace)./sqrt(num_cells2-1);
        shadedErrorBar_YS(plot_time, temp_trace_mean,temp_trace_sem, col_list(n_tr,:));
        axis tight;
        ylims1 = [min([ylims1(1); temp_trace_mean(:)-temp_trace_sem(:)]) max([ylims1(2); temp_trace_mean(:)+temp_trace_sem(:)])];
        if n_tr == 1
            ylabel(sprintf('%s, n=%d', ops.fig_title_run{n_flip}, num_cells2));
        else
            ax = gca;
            ax.YTickLabel = [];
            %ax.XTickLabel = [];
        end
        if n_flip < 3
            title(leg_list{n_pl});
        end
    end
end

for n_tr = 1:numel(plot_hand)
    ylim(plot_hand{n_tr}, ylims1);
end

end

