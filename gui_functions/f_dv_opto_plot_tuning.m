function f_dv_opto_plot_tuning(app)

data = app.data;

mice_all = unique(data.mouse_id, 'stable');

base_resp_win = [5 15];
params = f_dv_gather_params(app);

%%
num_dd_stim_dsets = 0;
peak_resp_all3 = cell(numel(mice_all),1);
mouse_id_all = cell(numel(mice_all),1);
for n_ms = 1:numel(mice_all)
    data2 = data(strcmpi(data.mouse_id, mice_all{n_ms}),:);
    fov_all = unique(data2.FOV_num, 'stable');
    
    peak_resp_all2 = cell(numel(fov_all),1);
    for n_fov = 1:numel(fov_all)
        data3 = data2(data2.FOV_num == fov_all(n_fov),:);
        
        idx_ammn_stim = strcmpi(data3.paradigm, 'ammn_stim');
        idx_ammn = strcmpi(data3.paradigm, 'ammn');
        peak_resp_all = cell(numel(data3.mouse_id),1);
        if sum(idx_ammn_stim)
            if sum(idx_ammn) > 1
                data_opto = data3(idx_ammn_stim, :);
                stim_ops = data_opto.proc_data{1}.stim_params.ops;
                if strcmpi(stim_ops.stim_trials_volt{2}, 'dev')
                    
                    cdata = cell(data_opto.num_planes,1);
                    for n_pl = 1:data_opto.num_planes
                        params.n_pl = n_pl;
                        cdata{n_pl} = f_dv_compute_cdata(data_opto, params);
                    end

                    cdata = cat(1, cdata{:});

                    firing_rate = cat(1,cdata.S_sm);

                    [~, T] = size(firing_rate);

                    chan_idx = strcmpi(data_opto.proc_ops{1}.chan_labels, 'Pockel');
                    stim_times = data_opto.proc_data{1}.stim_times_frame{chan_idx, n_pl};
                    stim_times(stim_times<base_resp_win(1)) = [];
                    stim_times(stim_times>(T-base_resp_win(2))) = [];

                    sorted_data = f_get_stim_trig_resp(firing_rate, stim_times, base_resp_win);

                    trial_ave = mean(sorted_data,3);
                    mean_resp = mean((trial_ave(:,(base_resp_win(1)+1):(base_resp_win(1)+8))),2);

        %             [~, idx2] = sort(mean_resp, 'descend');
        % 
        %             figure; imagesc(trial_ave(idx2,:));
        %             xlabel('frames'); ylabel('sorted cells');
        %             title('opto triggered response')
        %             
        %             figure; histogram(mean_resp)
        %             
                    thresh = mean(mean_resp) + 3*std(mean_resp);

                    stim_resp_cells = mean_resp>thresh;
                    
                    if sum(stim_resp_cells)
                        num_dd_stim_dsets = num_dd_stim_dsets + 1;
                    
                        stim_tn = stim_ops.stim_trials_volt{1} * 10;
                        stim_paradigm = stim_ops.stim_trials_volt{1};
                        mmn_freqs = stim_ops.MMN_patterns(stim_ops.paradigm_MMN_pattern(stim_paradigm),:);
                        stim_dd_freq = mmn_freqs(4-stim_paradigm);

                        cell_plane_lut = f_dv_opto_make_cell_plane_lut(data_opto);

                        resp_cells_dc_p_c = cell_plane_lut(stim_resp_cells,:);

                        resp_cells_dc_p_c_gc = [resp_cells_dc_p_c, zeros(size(resp_cells_dc_p_c,1),1)];
                        for n_cell = 1:size(resp_cells_dc_p_c_gc,1)
                            n_pl2 = resp_cells_dc_p_c(n_cell,2);
                            register_roi = [data_opto.register_roi{n_pl2}.fov_cell_idx, data_opto.register_roi{n_pl2}.reg_cell_idx];
                            n_cell2 = resp_cells_dc_p_c(n_cell,3);
                            idx4 = register_roi(:,2) == n_cell2;
                            resp_cells_dc_p_c_gc(n_cell,4) = register_roi(idx4,1);
                        end

                        stats1 = cat(1,data_opto.stats{:});
                        peak_val_all = cat(1,stats1.peak_vals);
                        peak_resp_all{idx_ammn_stim} = peak_val_all(resp_cells_dc_p_c_gc(:,1),stim_tn);

                        for n_dset2 = 1:numel(data3.mouse_id)
                            if idx_ammn(n_dset2)
                                data4 = data3(n_dset2,:);
                                stats2 = data4.stats;
                                num_cells2 = size(resp_cells_dc_p_c_gc,1);
                                peak_resp2 = nan(num_cells2,1);
                                for n_cell = 1:num_cells2
                                    n_pl2 = resp_cells_dc_p_c_gc(n_cell,2);
                                    register_roi = [data4.register_roi{n_pl2}.fov_cell_idx, data4.register_roi{n_pl2}.reg_cell_idx];
                                    n_gcell = resp_cells_dc_p_c_gc(n_cell,4);
                                    idx4 = register_roi(:,1) == n_gcell;
                                    n_cell2 = register_roi(idx4,2);
                                    if ~isnan(n_cell2)
                                        if 0%stats2{n_pl2}.peak_in_resp_win(n_cell2,stim_tn)
                                            peak_resp2(n_cell) = stats2{n_pl2}.peak_val_all(n_cell2,stim_tn);
                                        else
                                            peak_resp2(n_cell) = stats2{n_pl2}.trial_ave_lim_win_mean(n_cell2,1,stim_tn);
                                        end
                                    end
                                end
                                peak_resp_all{n_dset2} = peak_resp2;
                            end
                        end
                    end
                end
            end
        end
        peak_resp_all2{n_fov} = cat(2,peak_resp_all{:});
    end
    
    peak_resp_all3{n_ms} = cat(1,peak_resp_all2{:});
    mouse_id_all{n_ms} = ones(size(peak_resp_all3{n_ms},1),1)*n_ms;
end

peak_resp_all4 = cat(1, peak_resp_all3{:});
mouse_id_all2 = cat(1, mouse_id_all{:});

nan_free_idx = ~logical(sum(isnan(peak_resp_all4),2));
num_nanfree = sum(nan_free_idx);

nan_free_resp = peak_resp_all4(nan_free_idx,:);
nan_free_mouse = mouse_id_all2(nan_free_idx);
figure; hold on;
bar(mean(nan_free_resp), 'FaceColor', 'k', 'FaceAlpha', 0, 'LineWidth', 2)
errorbar(mean(nan_free_resp), std(nan_free_resp)/sqrt(size(nan_free_resp,1)-1), '.k', 'LineWidth', 2)
plot(nan_free_resp', 'o-', 'color', [.6 .6 .6]);
for n_ms = 1:numel(mice_all)
    idx1 = nan_free_mouse == n_ms;
    if sum(idx1)
        plot(mean(nan_free_resp(idx1,:),1)', 'o-', 'color', [0 0.4470 0.7410], 'LineWidth', 1)
    end
end
title('DD stimulated cell tuning'); ylabel('resp mag');

%% control dd ruun


num_dd_stim_dsets = 0;
peak_resp_all3 = cell(numel(mice_all),1);
mouse_id_all = cell(numel(mice_all),1);
for n_ms = 1:numel(mice_all)
    data2 = data(strcmpi(data.mouse_id, mice_all{n_ms}),:);
    fov_all = unique(data2.FOV_num, 'stable');
    
    peak_resp_all2 = cell(numel(fov_all),1);
    for n_fov = 1:numel(fov_all)
        data3 = data2(data2.FOV_num == fov_all(n_fov),:);
        
        idx_ammn_stim = strcmpi(data3.paradigm, 'ammn_stim');
        idx_ammn = strcmpi(data3.paradigm, 'ammn');
        peak_resp_all = cell(numel(data3.mouse_id),1);
        if sum(idx_ammn_stim)
            if sum(idx_ammn) > 1
                data_opto = data3(idx_ammn_stim, :);
                stim_ops = data_opto.proc_data{1}.stim_params.ops;
                if strcmpi(stim_ops.stim_trials_volt{2}, 'dev')
                    
                    cdata = cell(data_opto.num_planes,1);
                    for n_pl = 1:data_opto.num_planes
                        params.n_pl = n_pl;
                        cdata{n_pl} = f_dv_compute_cdata(data_opto, params);
                    end

                    cdata = cat(1, cdata{:});

                    firing_rate = cat(1,cdata.S_sm);

                    [~, T] = size(firing_rate);

                    chan_idx = strcmpi(data_opto.proc_ops{1}.chan_labels, 'Pockel');
                    stim_times = data_opto.proc_data{1}.stim_times_frame{chan_idx, n_pl};
                    stim_times(stim_times<base_resp_win(1)) = [];
                    stim_times(stim_times>(T-base_resp_win(2))) = [];

                    sorted_data = f_get_stim_trig_resp(firing_rate, stim_times, base_resp_win);

                    trial_ave = mean(sorted_data,3);
                    mean_resp = mean((trial_ave(:,(base_resp_win(1)+1):(base_resp_win(1)+8))),2);

        %             [~, idx2] = sort(mean_resp, 'descend');
        % 
        %             figure; imagesc(trial_ave(idx2,:));
        %             xlabel('frames'); ylabel('sorted cells');
        %             title('opto triggered response')
        %             
        %             figure; histogram(mean_resp)
        %             
                    thresh = mean(mean_resp) + 3*std(mean_resp);

                    stim_resp_cells = mean_resp>thresh;
                    
                    if sum(stim_resp_cells)
                        num_dd_stim_dsets = num_dd_stim_dsets + 1;
                    
                        stim_tn = stim_ops.stim_trials_volt{1} * 10;
                        stim_paradigm = stim_ops.stim_trials_volt{1};
                        mmn_freqs = stim_ops.MMN_patterns(stim_ops.paradigm_MMN_pattern(stim_paradigm),:);
                        stim_dd_freq = mmn_freqs(4-stim_paradigm);

                        cell_plane_lut = f_dv_opto_make_cell_plane_lut(data_opto);

                        resp_cells_dc_p_c = cell_plane_lut(stim_resp_cells,:);

                        resp_cells_dc_p_c_gc = [resp_cells_dc_p_c, zeros(size(resp_cells_dc_p_c,1),1)];
                        for n_cell = 1:size(resp_cells_dc_p_c_gc,1)
                            n_pl2 = resp_cells_dc_p_c(n_cell,2);
                            register_roi = [data_opto.register_roi{n_pl2}.fov_cell_idx, data_opto.register_roi{n_pl2}.reg_cell_idx];
                            n_cell2 = resp_cells_dc_p_c(n_cell,3);
                            idx4 = register_roi(:,2) == n_cell2;
                            resp_cells_dc_p_c_gc(n_cell,4) = register_roi(idx4,1);
                        end

                        stats1 = cat(1,data_opto.stats{:});
                        peak_val_all = cat(1,stats1.peak_vals);
                        peak_resp_all{idx_ammn_stim} = peak_val_all(resp_cells_dc_p_c_gc(:,1),stim_tn);

                        for n_dset2 = 1:numel(data3.mouse_id)
                            if idx_ammn(n_dset2)
                                data4 = data3(n_dset2,:);
                                stats2 = data4.stats;
                                num_cells2 = size(resp_cells_dc_p_c_gc,1);
                                peak_resp2 = nan(num_cells2,1);
                                for n_cell = 1:num_cells2
                                    n_pl2 = resp_cells_dc_p_c_gc(n_cell,2);
                                    register_roi = [data4.register_roi{n_pl2}.fov_cell_idx, data4.register_roi{n_pl2}.reg_cell_idx];
                                    n_gcell = resp_cells_dc_p_c_gc(n_cell,4);
                                    idx4 = register_roi(:,1) == n_gcell;
                                    n_cell2 = register_roi(idx4,2);
                                    if ~isnan(n_cell2)
                                        if 0%stats2{n_pl2}.peak_in_resp_win(n_cell2,stim_tn)
                                            peak_resp2(n_cell) = stats2{n_pl2}.peak_val_all(n_cell2,stim_tn);
                                        else
                                            peak_resp2(n_cell) = stats2{n_pl2}.trial_ave_lim_win_mean(n_cell2,1,stim_tn);
                                        end
                                    end
                                end
                                peak_resp_all{n_dset2} = peak_resp2;
                            end
                        end
                    end
                end
            end
        end
        peak_resp_all2{n_fov} = cat(2,peak_resp_all{:});
    end
    
    peak_resp_all3{n_ms} = cat(1,peak_resp_all2{:});
    mouse_id_all{n_ms} = ones(size(peak_resp_all3{n_ms},1),1)*n_ms;
end

peak_resp_all4 = cat(1, peak_resp_all3{:});
mouse_id_all2 = cat(1, mouse_id_all{:});

nan_free_idx = ~logical(sum(isnan(peak_resp_all4),2));
num_nanfree = sum(nan_free_idx);

nan_free_resp = peak_resp_all4(nan_free_idx,:);
nan_free_mouse = mouse_id_all2(nan_free_idx);
figure; hold on;
bar(mean(nan_free_resp), 'FaceColor', 'k', 'FaceAlpha', 0, 'LineWidth', 2)
errorbar(mean(nan_free_resp), std(nan_free_resp)/sqrt(size(nan_free_resp,1)-1), '.k', 'LineWidth', 2)
plot(nan_free_resp', 'o-', 'color', [.6 .6 .6]);
for n_ms = 1:numel(mice_all)
    idx1 = nan_free_mouse == n_ms;
    if sum(idx1)
        plot(mean(nan_free_resp(idx1,:),1)', 'o-', 'color', [0 0.4470 0.7410], 'LineWidth', 1)
    end
end
title('DD stimulated cell tuning'); ylabel('resp mag');


end