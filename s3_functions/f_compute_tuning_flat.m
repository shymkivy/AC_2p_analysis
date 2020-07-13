function data = f_compute_tuning_flat(data, ops)
    disp('Computing tuning...');
    
    for n_cond = 1:numel(ops.conditions_to_analyze)
        cond_name = ops.conditions{ops.conditions_to_analyze(n_cond)}; 
        num_cells = data.(cond_name).num_cells;
        num_frames = size(data.(cond_name).vec_data,2);

        trial_ave = zeros(num_cells, num_frames, numel(ops.context_types_all));
        
        % for all cells find averages over trials and subtract baselines
        for n_cell = 1:num_cells
            % there are two sources of variance, the gaussian white noise, and
            % the signal variance. when whe average over trials, the white
            % noise averages out, but the signal should remain. when we average
            % over trials, we shoud not divide by std, because the white noise
            % is different for every cell. To normalize everything, we adjusted
            % that peak firing rate to 1 previously.

            % pull out data
            if ops.use_spike_inf
                data_load = squeeze(data.(cond_name).vec_data_si(n_cell,:,:));
            else
                data_load = squeeze(data.(cond_name).vec_data(n_cell,:,:));
            end
            trial_type = data.(cond_name).trial_type(n_cell,:);
            
            
            % average
            for n_cntxt = 1:numel(ops.context_types_all)
                if sum(trial_type == ops.context_types_all(n_cntxt)) > 0
                    trial_ave(n_cell,:,n_cntxt) = mean(data_load(:,trial_type == ops.context_types_all(n_cntxt)),2);
                else
                    trial_ave(n_cell,:,n_cntxt) = zeros(1,num_frames);
                end
            end

            % find the baseline to subtract. Use KDE and find peak
            [f,xi] = ksdensity(reshape(trial_ave(n_cell,:,:),[],1));
            baseline = mean(xi(f == max(f)));  % added mean not sure if is ok
            trial_ave(n_cell,:,:) = trial_ave(n_cell,:,:) - baseline;

            % plot KDE
    %         figure
    %         plot(xi,f);
    %         title(sprintf('cell %d',n_cell))

    %         % plots for fun
    %         figure;
    %         for ii = 1:numel(ops.context_types_all)
    %             subplot(4,7,ii)
    %             %plot(ops.time_stim_window, temp_data(:,temp_trial_type == ops.context_types_all(ii)))
    %             hold on;
    %             plot(ops.time_stim_window, trial_ave(n_cell,:,ii) , 'k', 'LineWidth', 2);
    %             ylim([-0.1 0.3]);
    %             title(context_type2(ii));
    %         end
    %         suptitle(sprintf('cell %d',n_cell))
        end

        z_factor = std(reshape(trial_ave, [], 1));

    %     % abnother way
    %     [f,xi] = ksdensity(reshape(trial_ave(:,:,:),[],1));
    %     z_factor = rms(trial_ave(trial_ave<xi(f == max(f))));
        
        trial_ave_z = trial_ave/z_factor;

        % find cells/trials that cross z threshold at each context
        resp_cells = false(num_cells, numel(ops.context_types_all));
        resp_cells_zmag = zeros(num_cells, numel(ops.context_types_all));
        for n_cell = 1:num_cells
            for n_cntxt = 1:numel(ops.context_types_all)
                resp_cells(n_cell, n_cntxt) = sum(trial_ave_z(n_cell, ops.win.resp_window_frames,n_cntxt) > ops.z_threshold)>0;
                resp_cells_zmag(n_cell, n_cntxt) = mean(trial_ave_z(n_cell, ops.win.resp_window_frames,n_cntxt));
            end
        end

        % find all tuned cells
        tuned_cells = sum(resp_cells(:,1:10),2)>0;
        % how many cells are tuned to each freq


    %     % plots for fun of responses adn thresholds
    %     y_max = max(trial_ave_z(:));
    %     y_min = min(trial_ave_z(:));
    %     for n_cell = 1:15
    %         if tuned_cells(n_cell)
    %             figure;
    %             for ii = 1:numel(ops.context_types_all)
    %                 subplot(4,7,ii)
    %                 %plot(ops.time_stim_window, temp_data(:,temp_trial_type == ops.context_types_all(ii)))
    %                 hold on;
    %                 plot(ops.time_stim_window, trial_ave_z(n_cell,:,ii) , 'k', 'LineWidth', 2);
    %                 l1 = line([ops.time_stim_window(1) ops.time_stim_window(end)],[z_threshold z_threshold], 'Color', 'r');
    %                 ylim([y_min y_max]);
    %                 title(ops.context_types_all(ii));
    %             end
    %             suptitle(sprintf('cell %d',n_cell))
    %         end
    %     end




        % create pattern key for processing for every cell (control index, red index, dev index, ...)
        context_mmn = zeros(num_cells, 3, 2);
        for n_cell = 1:num_cells
            temp_MMN_freq = data.(cond_name).MMN_freq(n_cell,:);

            % use ops.context_types_all as key
            %context_mmn(n_cell,:,:) = [temp_MMN_freq(2), 201 + params.redundent_to_analyze, 170; ...
            %                           temp_MMN_freq(1), 100 + params.redundent_to_analyze, 270]';
            context_mmn(n_cell,:,:) = [temp_MMN_freq(2), 10 + ops.redundent_to_analyze, 19; ...
                                       temp_MMN_freq(1), 19 + ops.redundent_to_analyze, 28]';
        end

        % save responses for the mmn trials for all cells and find responsive cells
        trial_ave_z_mmn = zeros(num_cells, num_frames, 3, 2);
        trial_ave_mmn = zeros(num_cells, num_frames, 3, 2);
        context_resp_cells = false(num_cells, 2);
        for n_cell = 1:num_cells
            for n_flip = 1:2
                % pull out traces of contextual trials (cont, red, dev)
                trial_ave_z_mmn(n_cell, :, :, n_flip) = trial_ave_z(n_cell, : ,context_mmn(n_cell,:,n_flip));
                trial_ave_mmn(n_cell, :, :, n_flip) = trial_ave(n_cell, : ,context_mmn(n_cell,:,n_flip));
                % mark responsive cell if cell responds in any of the 3 contexts
                context_resp_cells(n_cell, n_flip) = sum(reshape(trial_ave_z(n_cell, ops.win.resp_window_frames,context_mmn(n_cell,:,n_flip)), [],1) > ops.z_threshold)>0;
            end
        end
%         figure; plot(mean(trial_ave_z_mmn(:,:,1,1),1)); hold on;
%         plot(mean(trial_ave_z_mmn(:,:,2,1),1));
%         plot(mean(trial_ave_z_mmn(:,:,3,1),1));
% 
%         figure; plot(mean(trial_ave_z_mmn(:,:,1,2),1)); hold on;
%         plot(mean(trial_ave_z_mmn(:,:,2,2),1));
%         plot(mean(trial_ave_z_mmn(:,:,3,2),1));
%         
%         
%         figure; hold on;
%         plot(mean(trial_ave_z(:,:,7),1));
%         plot(mean(trial_ave_z(:,:,22),1));
%         plot(mean(trial_ave_z(:,:,19),1));
%         
%         figure; hold on;
%         plot(mean(trial_ave_z(:,:,5),1));
%         plot(mean(trial_ave_z(:,:,13),1));
%         plot(mean(trial_ave_z(:,:,28),1));
        

        % number of cells responding in each context
        %sum(context_resp_cells)

        % also compute axis limits
        
        y_min = zeros(2,1);
        y_max = zeros(2,1);
        for n_flip = 1:2
            temp_mean = reshape(mean(trial_ave_z_mmn(context_resp_cells(:,n_flip),:,:,n_flip)),[],1);
            temp_sem = reshape(std(trial_ave_z_mmn(context_resp_cells(:,n_flip),:,:,n_flip)),[],1)/sqrt(sum(context_resp_cells(:,n_flip))-1);
            y_min(n_flip) = min(temp_mean - temp_sem);
            y_max(n_flip) = max(temp_mean + temp_sem);
        end
        
        %figure; imagesc(trial_ave(:,:,4))
        
        % plot psth
%         figure; plot(ops.time_stim_window(5:end), mean(trial_ave(:,5:end,9),1));

%         figure; hold on;
%         plot(ops.time_stim_window(5:end), mean(trial_ave_z(:,5:end,7),1), 'k');
%         plot(ops.time_stim_window(5:end), mean(trial_ave_z(:,5:end,22),1), 'b');
%         plot(ops.time_stim_window(5:end), mean(trial_ave_z(:,5:end,19),1), 'r');
%         
%         figure; hold on;
%         plot(ops.time_stim_window(5:end), mean(trial_ave_z(:,5:end,5),1), 'k');
%         plot(ops.time_stim_window(5:end), mean(trial_ave_z(:,5:end,13),1), 'b');
%         plot(ops.time_stim_window(5:end), mean(trial_ave_z(:,5:end,28),1), 'r');
%         

        % save
        data.(cond_name).trial_ave_z = trial_ave_z;
        data.(cond_name).z_factor = z_factor;
        data.(cond_name).z_threshold = ops.z_threshold;
        data.(cond_name).resp_cells = resp_cells;
        data.(cond_name).tuned_cells = tuned_cells;
        data.(cond_name).tuned_cells_totals = sum(resp_cells(tuned_cells,1:10));
        data.(cond_name).trial_ave_mmn = trial_ave_mmn;
        data.(cond_name).trial_ave_z_mmn = trial_ave_z_mmn;
        data.(cond_name).context_resp_cells = context_resp_cells;
        data.(cond_name).num_context_resp_cells = sum(context_resp_cells);
        data.(cond_name).context_mmn = context_mmn;
        data.(cond_name).y_min = y_min;
        data.(cond_name).y_max = y_max;
        data.(cond_name).resp_cells_zmag = resp_cells_zmag;
    end
    
    y_min = zeros(numel(ops.conditions_to_analyze),2);
    y_max = zeros(numel(ops.conditions_to_analyze),2);
    for n_cond = 1:numel(ops.conditions_to_analyze)
        cond_name = ops.conditions{ops.conditions_to_analyze(n_cond)}; 
        y_min(n_cond,:) = data.(cond_name).y_min;
        y_max(n_cond,:) = data.(cond_name).y_max;
    end
    data.gen_params.y_min = min(y_min(:));
    data.gen_params.y_max = max(y_max(:));
end