function f_mpl_plot_ctx_cond(data, ops)

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    
    resp_cells = logical(sum(cat(1,cdata.peak_tuned_trials_combined_ctx{:}),2));
    trial_ave = cat(1,cdata.trial_ave{:});
 
    MMN_freq_all = cell(numel(cdata.area),1);
    for n_dset = 1:numel(cdata.area)
        MMN_freq_all{n_dset} = repmat(cdata.MMN_freq{n_dset},cdata.num_cells(n_dset),1);
    end
    MMN_freq_all = cat(1,MMN_freq_all{:});
   
    plot_traces = trial_ave(resp_cells,:,:);
    f_mpl_plot_ctx(plot_traces, MMN_freq_all(resp_cells,:), cdata.trial_window{n_dset}.trial_window_t, ops);
    suptitle(sprintf('%s, %d resp cells', cond_name, size(plot_traces,1)));
        
    
%     for n_dset = 1:data.(cond_name).num_dsets
%         trial_ave_z = data.(cond_name).trial_ave_z{n_dset};
%         tuned_cells = logical(sum(data.(cond_name).resp_cells{n_dset}(:,1:10),2));
%         context_resp_cells = data.(cond_name).resp_cells_mmn{n_dset};
% 
%         MMN_freq = data.(cond_name).MMN_freq{n_dset};
%         trial_window_t = data.(cond_name).trial_window_t{n_dset};
%         plot_time = trial_window_t;
%         
%         plot_cells = trial_ave_z(:,:,:);
%         plot_trace = squeeze(mean(plot_cells,1)); 
%         plot_contx_traces(plot_time, plot_trace, MMN_freq);
%         suptitle(sprintf('%s, all cells dataset %d, %d cells', cond_name, n_dset, size(plot_cells,1)));
%         
%         plot_cells = trial_ave_z(tuned_cells,:,:);
%         plot_trace = squeeze(mean(plot_cells,1));
%         plot_contx_traces(plot_time, plot_trace, MMN_freq)
%         suptitle(sprintf('%s, freq tuned cells %d, %d cells', cond_name, n_dset,size(plot_cells,1)));
%         
%         plot_cells = trial_ave_z(context_resp_cells(:,1),:,:);
%         plot_trace = squeeze(mean(plot_cells,1));
%         plot_contx_traces(plot_time, plot_trace, MMN_freq)
%         suptitle(sprintf('%s, mmn tuned cells dataset %d, %d cells', cond_name, n_dset,size(plot_cells,1)));
%         
%         plot_cells = trial_ave_z(context_resp_cells(:,2),:,:);
%         plot_trace = squeeze(mean(plot_cells,1));
%         plot_contx_traces(plot_time, plot_trace, MMN_freq)
%         suptitle(sprintf('%s, flipmmn tuned cells dataset %d, %d cells', cond_name, n_dset,size(plot_cells,1)));
%     end
end

end

% function plot_contx_traces(plot_time, plot_trace, MMN_freq)
% ylims1 = [min(plot_trace(:)) max(plot_trace(:))];
% figure;
% 
% % first line
% for ii = 1:10
%     subplot(3,10,ii)
%     cute_plot(plot_time, plot_trace(:,ii));
%     ylim(ylims1);
%     title(sprintf('freq %d', ii));
%     ax = gca;
%     if ii == 1
%         temp_tick = ax.XTick;
%     elseif ii > 1
% %         ax.YTickLabel = [];
% %         ax.XTickLabel = [];
%     end
% end
% % second line
% subplot(3,10,11);
% cute_plot(plot_time, plot_trace(:,MMN_freq(2)));
% ylim(ylims1);
% title(sprintf('cont'));
% ylabel(sprintf('freq %d', MMN_freq(2)))
% ax = gca;
% ax.XTick = temp_tick;
% for ii = 11:18
%     subplot(3,10,ii+1);
%     cute_plot(plot_time, plot_trace(:,ii));
%     ylim(ylims1);
%     title(sprintf('red %d', ii-10));
%     ax = gca;
% %     ax.YTickLabel = [];
% %     ax.XTickLabel = [];
% end
% subplot(3,10,20);
% cute_plot(plot_time, plot_trace(:,19));
% ylim(ylims1);
% title(sprintf('dev'));
% ax = gca;
% % ax.YTickLabel = [];
% % ax.XTickLabel = [];
% % line 3
% subplot(3,10,21);
% cute_plot(plot_time, plot_trace(:,MMN_freq(1)));
% ylim(ylims1);
% title(sprintf('cont'));
% ylabel(sprintf('freq %d', MMN_freq(1)))
% ax = gca;
% ax.XTick = temp_tick;
% for ii = 20:27
%     subplot(3,10,ii+2);
%     cute_plot(plot_time, plot_trace(:,ii));
%     ylim(ylims1);
%     title(sprintf('red %d', ii-19));
%     ax = gca;
% %     ax.YTickLabel = [];
% %     ax.XTickLabel = [];
% end
% subplot(3,10,30);
% cute_plot(plot_time, plot_trace(:,28));
% ylim(ylims1);
% title(sprintf('dev'));
% ax = gca;
% % ax.YTickLabel = [];
% % ax.XTickLabel = [];
% 
% end
% 
% function cute_plot(plot_time, plot_trace)
%     plot(plot_time, plot_trace, 'Color', [1 0 1], 'LineWidth', 2);
% end
