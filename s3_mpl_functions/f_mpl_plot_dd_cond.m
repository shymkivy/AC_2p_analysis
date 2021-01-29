function f_mpl_plot_dd_cond(data, ops)
%     y_min = data.gen_params.y_min;
%     y_max = data.gen_params.y_max;

sp = cell(numel(ops.regions_to_analyze)*numel(ops.flip_to_analyze),1);
dev_ymin = 0;
dev_ymax = 0;
figure;
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);

    resp_cells_mmn = logical(cat(1,cdata.peak_tuned_trials_combined_ctx{:}));
    trial_ave_mmn = cat(1,cdata.trial_ave_mmn{:});
    trial_window_t = cdata.trial_window{1}.trial_window_t;
    
    plot_traces = cell(3,1);
    if strcmpi(ops.dev_cells_ctx, 'ctx_tuned')
        plot_cells = logical(sum(resp_cells_mmn(:,1:3),2));
        plot_traces{1} = trial_ave_mmn(plot_cells,:,1:3);
        plot_cells = logical(sum(resp_cells_mmn(:,4:6),2));
        plot_traces{2} = trial_ave_mmn(plot_cells,:,4:6);
    elseif strcmpi(ops.dev_cells_ctx, 'tuned_all')
        plot_cells = logical(sum(cat(1,cdata.resp_cells_all_onset{:}) + cat(1,cdata.resp_cells_all_offset{:}),2));
        plot_traces{1} = trial_ave_mmn(plot_cells,:,1:3);
        plot_traces{2} = trial_ave_mmn(plot_cells,:,4:6);
    elseif strcmpi(ops.dev_cells_ctx, 'all')
        plot_traces{1} = trial_ave_mmn(:,:,1:3);
        plot_traces{2} = trial_ave_mmn(:,:,4:6);
    end
    plot_traces{3} = cat(1,plot_traces{1:2});
    
    for n_flip = 1:3
  
        
%         % plot individual cells
%         close all
%         for n_cell = 1:30
%             if temp_resp_cells(n_cell)
%                 figure;
%                 hold on;
%                 plot(ops.time_stim_window, temp_data(n_cell,:,1), 'k');
%                 plot(ops.time_stim_window, temp_data(n_cell,:,2), 'b');
%                 plot(ops.time_stim_window, temp_data(n_cell,:,3), 'r');
%                 l1 = line([ops.time_stim_window(1) ops.time_stim_window(end)],[params.z_threshold params.z_threshold], 'Color', 'r');
%                 title(sprintf('Cell %d', n_cell));
%             end
%         end

        traces_ave = squeeze(mean(plot_traces{n_flip}));
        traces_SEM = squeeze(std(plot_traces{n_flip}))/sqrt(size(plot_traces{n_flip},1)-1);
        dev_ymin = min([dev_ymin; traces_ave(:)-traces_SEM(:)]);
        dev_ymax = max([dev_ymax; traces_ave(:)+traces_SEM(:)]);
%             figure; hold on;
%             plot(traces_ave(:,1))
%             plot(traces_ave(:,2))
%             plot(traces_ave(:,3))

        % plot limited trials data
        plot_ind = n_cond+(n_flip-1)*numel(ops.regions_to_analyze);
        sp{plot_ind} = subplot(numel(ops.flip_to_analyze),numel(ops.regions_to_analyze),plot_ind);
        %figure;
        hold on;
        if size(plot_traces{n_flip},1)
            for n_cntxt = 1:3
                % converting to z-scores is tricky, I dont do it here yet
                %patch([0 0 ops.stim_duration ops.stim_duration],[y_min y_max y_max y_min] ,[224, 243, 255]/256, 'LineStyle', 'none', 'FaceAlpha', 0.5);
                shadedErrorBar_YS(trial_window_t, traces_ave(:,n_cntxt), traces_SEM(:,n_cntxt), ops.context_colors{n_cntxt});

            end
        end
        axis tight;
        if n_cond == 1
            ylabel('z-score');
            if n_flip == 3
                xlabel('time (sec)');
            end
        end
        title(sprintf('%s, %s, %d cells', cond_name, ops.fig_title_run{n_flip}, size(plot_traces{n_flip},1)));   
    end 
end
suptitle(sprintf('%s: z-thresh=%d; redundant num %d', ops.paradigm_type, ops.stat.z_scores_thresh, ops.redundent_to_analyze));
for n_sp = 1:numel(sp)
    ylim(sp{n_sp}, [dev_ymin dev_ymax]);
end


end