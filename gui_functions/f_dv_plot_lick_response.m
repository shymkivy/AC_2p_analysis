function f_dv_plot_lick_response(app)

ddata = app.ddata;

if ~isempty(ddata.stats_within{1})
    trial_data_mean_z = cell(ddata.num_planes,1);
    trial_data_mean_z_rew = cell(ddata.num_planes,1);
    trial_data_mean_z_notrew = cell(ddata.num_planes,1);
    
    stat_window_t = ddata.stats_within{1}.stat_window_t;
    win_use = [-0.5 8];
    win_on = [0 0.5];
    win_off = [0.5 2];
    idx_win_use = logical((stat_window_t > win_use(1)).* (stat_window_t < win_use(2)));
    
    onset_idx = logical((stat_window_t > win_on(1)).* (stat_window_t < win_on(2)));
    offset_idx = logical((stat_window_t > win_off(1)).* (stat_window_t < win_off(2)));
    
    for n_pl = 1:ddata.num_planes
        stats1 = ddata.stats_within{n_pl};
        
        trial_data_mean_z{n_pl} = stats1.stim_trial_data_mean_z;
        trial_data_mean_z_rew{n_pl} = stats1.trial_data_mean_z_rew;
        trial_data_mean_z_notrew{n_pl} = stats1.trial_data_mean_z_notrew;
    end
    
    trial_data_mean_z = cat(1, trial_data_mean_z{:});
    trial_data_mean_z_rew = cat(1, trial_data_mean_z_rew{:});
    trial_data_mean_z_notrew = cat(1, trial_data_mean_z_notrew{:});
    
    [peak_val, peak_loc] = max(trial_data_mean_z, [], 2);
    [peak_val_on, peak_loc_on] = max(trial_data_mean_z(:, onset_idx), [], 2);
    [peak_val_off, peak_loc_off] = max(trial_data_mean_z(:, offset_idx), [], 2);
    [peak_val_rew, peak_loc_rew] = max(trial_data_mean_z_rew, [], 2);
    [peak_val_notrew, peak_loc_notrew] = max(trial_data_mean_z_notrew, [], 2);
    
    [~, cell_ord] = sort(peak_val, 'descend');
    [~, cell_ord_on] = sort(peak_val_on, 'descend');
    [~, cell_ord_off] = sort(peak_val_off, 'descend');
    [~, cell_ord_rew] = sort(peak_val_rew, 'descend');
    [~, cell_ord_notrew] = sort(peak_val_notrew, 'descend');
    
    peak_vals_all = [peak_val, peak_val_on, peak_val_off, peak_val_rew, peak_val_notrew];
    %peak_vals_all = [peak_val_off, peak_val_rew, peak_val_notrew];
    %peak_vals_all(peak_vals_all<6) = 0;
    
    % sort nonrewarded
%     [~, cell_ord1] = sort(peak_vals_all(:,1), 'descend');
%     [~, cell_ord2] = sort(peak_vals_all(:,2), 'descend');
%     [~, cell_ord3] = sort(peak_vals_all(:,3), 'descend');
%     [~, cell_ord4] = sort(peak_vals_all(:,4), 'descend');
%     [~, cell_ord5] = sort(peak_vals_all(:,5), 'descend');
% 
%     [~, cell_ord6] = sort(max(peak_vals_all, [], 2), 'descend');
    
    keep_idx = logical(sum(peak_vals_all>6,2));
    peak_vals_alln = peak_vals_all(keep_idx,:);
    peak_vals_alln = peak_vals_alln - mean(peak_vals_alln);
    peak_vals_alln = peak_vals_alln./std(peak_vals_alln);
    
    params = struct();
    hc_out = f_hcluster_wrap(peak_vals_alln, params);
    
    ord_all = hc_out.dend_order;
    
    trial_data_mean_z2 = trial_data_mean_z(keep_idx,:);
    trial_data_mean_z_rew2 = trial_data_mean_z_rew(keep_idx,:);
    trial_data_mean_z_notrew2 = trial_data_mean_z_notrew(keep_idx,:);
    
    clim1 = [-6 10];
    figure;
    subplot(1, 3, 1);
    imagesc(stat_window_t(idx_win_use), [], trial_data_mean_z2(ord_all,idx_win_use))
    caxis(clim1)
    title('tone responsive');
    ylabel('Cells z-scored, sorted');
    xlabel('Time (sec)');
    subplot(1, 3, 2);
    imagesc(stat_window_t(idx_win_use), [], trial_data_mean_z_rew2(ord_all,idx_win_use))
    title('lick - reward');
    caxis(clim1)
    subplot(1, 3, 3);
    imagesc(stat_window_t(idx_win_use), [], trial_data_mean_z_notrew2(ord_all,idx_win_use))
    title('no lick no rewarded');
    caxis(clim1)
%     
%     figure; 
%     imagesc(stat_window_t(idx_win_use), [], trial_data_mean_z(cell_ord_on,idx_win_use))
%     caxis([-6 10])
%     
%     figure; 
%     imagesc(stat_window_t(idx_win_use), [], trial_data_mean_z(cell_ord_off(),idx_win_use))
%     caxis([-6 10])
%     
%     figure; 
%     imagesc(stat_window_t(idx_win_use), [], trial_data_mean_z_rew(cell_ord_rew,idx_win_use))
%     caxis([-6 10])
%     
%     figure; 
%     imagesc(stat_window_t(idx_win_use), [], trial_data_mean_z_notrew(cell_ord_notrew,idx_win_use))
%     caxis([-6 10])
    

    stim_times_volt = ddata.proc_data{1}.stim_times_volt{1}/1000;
    lick_times_volt = ddata.proc_data{1}.stim_times_volt{5}/1000;
    rew_times_volt = ddata.proc_data{1}.stim_times_volt{6}/1000;
    
    win_size = [-1, 5];
    
    win_t = win_size(1):0.01:win_size(2);
    
    num_stim = numel(stim_times_volt);
    lick_raster = zeros(num_stim, numel(win_t));
    
    for n_st = 1:num_stim
        st_time = stim_times_volt(n_st);
        temp1 = lick_times_volt - st_time;
        idx1 = logical((temp1 > win_size(1)) .* (temp1 < win_size(2)));
        licks1 = temp1(idx1);
        num_licks2 = numel(licks1);
        for n_lick = 1:num_licks2
            idx2 = licks1(n_lick) < win_t;
            if sum(idx2)
                idx3 = find(idx2,1);
                lick_raster(n_st, idx3) = 1;
            end
        end

    end
    
    figure; hold on; 
    im1 = imagesc(win_t, [], 1-lick_raster);
    line([0 0], [0 num_stim+1], 'color', 'r')
    line([.5 .5], [0 num_stim+1], 'color', 'r')
    colormap('gray'); axis tight;
    im1.Parent.YDir =  'Reverse';
    title('Lick')
    ylabel('Trial');
    xlabel('Time (sec)');
    
end
    
%[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

idx1 = logical(strcmpi(ddata.mouse_id, app.data.mouse_id).*(ddata.FOV_num == app.data.FOV_num));
data3 = app.data(idx1,:);


end