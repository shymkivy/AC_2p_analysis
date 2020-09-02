function f_mpl_plot_ctx4(trial_ave, resp_cells,plot_time,ops)
%% this will plot 1-n redundent and the dev
if ops.plot_combined
    m_plt = 3;
else
    m_plt = 2;
end

%%
plt_sequence = [8, 1:7, 10];
num_n_traces = numel(plt_sequence);

%%
col_list = [ops.context_colors{1}; repmat(ops.context_colors{2}, 7, 1); ops.context_colors{3}];

resp_cells1 = reshape(resp_cells, [], 3,2);

%leg_list = [ops.context_types_labels(plt_sequence+10); ops.context_types_labels(plt_sequence+20)];   

leg_list = {'Cont', 'RedF Dev', 'ContF', 'Red DevF'};   

%% get traces
traces1{1} = trial_ave(:,:,plt_sequence+10);
traces1{2} = trial_ave(:,:,plt_sequence+20);

traces1{1} = traces1{1}(logical(sum(resp_cells1(:,:,1),2)),:,:);
traces1{2} = traces1{2}(logical(sum(resp_cells1(:,:,2),2)),:,:);

if ops.plot_combined
    traces1{3} = cat(1, traces1{1}, traces1{2});
end

%% plot

plot_hand = cell(m_plt*2,1);
ylims1 = [0 0.2];

figure;
for n_flip = 1:m_plt
    traces2 = traces1{n_flip};
    num_cells2 = size(traces2,1);
    time_shift = 0;
    for n_tr = 1:num_n_traces
        n_pl_shift = num_n_traces*(n_flip-1);
        if n_tr == 1
            plot_hand{2*(n_flip-1)+n_tr} = subplot(m_plt,num_n_traces,1+n_pl_shift);
        elseif n_tr == 2
            plot_hand{2*(n_flip-1)+n_tr} = subplot(m_plt,num_n_traces,(2:num_n_traces)+n_pl_shift); hold on;
        end
        temp_trace = traces2(:,:,n_tr);
        temp_trace_mean = mean(temp_trace);
        temp_trace_sem = std(temp_trace)./sqrt(num_cells2-1);
        shadedErrorBar_YS(plot_time+time_shift, temp_trace_mean,temp_trace_sem, col_list(n_tr,:));
        axis tight;
        ylims1 = [min([ylims1(1); temp_trace_mean(:)-temp_trace_sem(:)]) max([ylims1(2); temp_trace_mean(:)+temp_trace_sem(:)])];
        if n_tr == 1
            ylabel(sprintf('%s, n=%d', ops.fig_title_run{n_flip}, num_cells2));
        else
            time_shift = time_shift+1;
            ax = gca;
            ax.YTickLabel = [];
            %ax.XTickLabel = [];
        end
        if n_flip < 3
            if n_tr == 1 || n_tr ==2
                title(leg_list{2*(n_flip-1)+n_tr});
            end
        end
    end
end

for n_tr = 1:numel(plot_hand)
    ylim(plot_hand{n_tr}, ylims1);
end

end