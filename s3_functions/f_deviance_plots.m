function f_deviance_plots(data, ops)
    y_min = data.gen_params.y_min;
    y_max = data.gen_params.y_max;
    figure;
    for n_flip = 1:3
        for n_cond = 1:numel(ops.conditions_to_analyze)
            cond_name = ops.conditions{ops.conditions_to_analyze(n_cond)};

            if n_flip == 3
                % combine flip
                ctx_cells = [data.(cond_name).context_resp_cells(:,1);...
                                   data.(cond_name).context_resp_cells(:,2)];
                               
                if ~ops.dev_ctx_cells
                    ctx_cells = ones(size(ctx_cells,1),size(ctx_cells,2));
                end
                
                if ops.dev_use_z
                    ctx_traces = [data.(cond_name).trial_ave_z_mmn(:,:,:,1);...
                                data.(cond_name).trial_ave_z_mmn(:,:,:,2)];
                else
                    ctx_traces = [data.(cond_name).trial_ave_mmn(:,:,:,1);...
                                data.(cond_name).trial_ave_mmn(:,:,:,2)];
                end
                ctx_traces = ctx_traces(ctx_cells,:,:);
            else

                %(num_cells, frames, context, flip)
                ctx_cells = data.(cond_name).context_resp_cells(:,n_flip);
                if ~ops.dev_ctx_cells
                    ctx_cells = ones(size(ctx_cells,1),size(ctx_cells,2));
                end
                
                if ops.dev_use_z
                    ctx_traces = data.(cond_name).trial_ave_z_mmn(ctx_cells,:,:,n_flip);
                else
                    ctx_traces = data.(cond_name).trial_ave_mmn(ctx_cells,:,:,n_flip);
                end
            end
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

            traces_ave = squeeze(mean(ctx_traces));
            traces_SEM = squeeze(std(ctx_traces, [], 1))/sqrt(sum(ctx_cells)-1);
            
%             figure; hold on;
%             plot(traces_ave(:,1))
%             plot(traces_ave(:,2))
%             plot(traces_ave(:,3))
            
            % plot limited trials data
            subplot(3,4,n_cond+(n_flip-1)*4)
            %figure;
            hold on;
            for n_cntxt = 1:3
                % converting to z-scores is tricky, I dont do it here yet
                %patch([0 0 ops.stim_duration ops.stim_duration],[y_min y_max y_max y_min] ,[224, 243, 255]/256, 'LineStyle', 'none', 'FaceAlpha', 0.5);
                shadedErrorBar(ops.time_stim_window(ops.win.analysis_window), (traces_ave(ops.win.analysis_window,n_cntxt)), traces_SEM(ops.win.analysis_window,n_cntxt), 'lineprops', ops.context_colors{n_cntxt});

            end
            axis tight;
            ylim([y_min y_max])
            if n_cond == 1
                ylabel('z-score');
                if n_flip == 3
                    xlabel('time (sec)');
                end
            end
            title(sprintf('%s, %s', ops.conditions{n_cond}, ops.fig_title_run{n_flip}));   
        end 
    end
    suptitle(sprintf('%s: z-thresh=%d; redundant num %d', ops.paradigm_type, ops.z_threshold, ops.redundent_to_analyze));

end