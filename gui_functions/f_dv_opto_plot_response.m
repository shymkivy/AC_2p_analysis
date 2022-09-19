function f_dv_opto_plot_response(app)

n_pl = app.mplSpinner.Value;

ddata = app.ddata;

cdata = f_dv_get_cdata(app);
firing_rate = cat(1,cdata.S_sm);

cdata.volume_period


trial_window = [-.5, 2];
stim_analysis_win = [.1 .4];

chan_idx = strcmpi(ddata.proc_ops{1}.chan_labels, 'Pockel');
stim_times = ddata.proc_data{1}.stim_times_frame{chan_idx, n_pl};

[trial_window_t, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1,:).volume_period);
stim_analysis_idx = and(trial_window_t >= stim_analysis_win(1), trial_window_t <= stim_analysis_win(2));

fprintf('%d laser stim, every %.2f sec\n', numel(stim_times), mean(diff(stim_times)))

if ~isempty(stim_times)

    [~, T] = size(firing_rate);

    stim_times(stim_times<trial_frames(1)) = [];
    stim_times(stim_times>(T-trial_frames(2))) = [];

    sorted_data = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
    
    % for inserting color bar under imagesc
    %opto_trace = zeros(1, sum(trial_frames), 3);
    %opto_trace(1,trial_window_t==0,1) = 1;

    trial_ave = mean(sorted_data,3);
    mean_resp = mean(trial_ave(:,stim_analysis_idx),2);

    [mean_resp_sort, idx1] = sort(mean_resp, 'descend');

    f1 = figure;hold on; 
    imagesc(trial_window_t, [], trial_ave(idx1,:));
    line([0 0], [size(trial_ave,1)+0.5, 0.5], 'color', [1 0 0], 'LineWidth', 1)
    %imagesc(trial_window_t, size(trial_ave,1)+1, opto_trace);
    axis tight;
    f1.Children.YDir = 'reverse';
    xlabel('Time (sec)'); ylabel('Sorted cells');
    title('Opto triggered response')
    
    sort_data2 = [sorted_data(idx1(1),:,:), sorted_data(idx1(end),:,:)];
    ylims1 = [min(sort_data2(:)) max(sort_data2(:))];
    
    figure; hold on;
    plot(trial_window_t, squeeze(sorted_data(idx1(1),:,:)), 'color', [.6 .6 .6]);
    plot(trial_window_t, mean(sorted_data(idx1(1),:,:),3), 'm', 'linewidth', 2);
    line([0 0], ylims1, 'color', [1 0 0], 'LineWidth', 1)
    %imagesc(trial_window_t, -1, opto_trace);
    axis tight; xlabel('Time (sec)'); ylabel('response magnitude');
    title('Top activated cell')
    ylim(ylims1); axis tight
    
    figure; hold on;
    plot(trial_window_t, squeeze(sorted_data(idx1(end),:,:)), 'color', [.6 .6 .6]);
    plot(trial_window_t, mean(sorted_data(idx1(end),:,:),3), 'm', 'linewidth', 2);
    line([0 0], ylims1, 'color', [1 0 0], 'LineWidth', 1)
    axis tight; xlabel('Time (sec)'); ylabel('response magnitude');
    title('Top inhibited cell')
    ylim(ylims1);  axis tight
    
    [num_cells, ~, num_trials] = size(sorted_data);
    mean_resp_by_trial = squeeze(mean(sorted_data(idx1, stim_analysis_idx, :),2));
    
    x_data = 1:num_trials;
    
    y1_all = zeros(num_cells,1);
    y0_all = zeros(num_cells,1);
    fitmeans = zeros(num_cells,1);
    for n_cell = 1:num_cells
        Y_data = mean_resp_by_trial(n_cell,:);
        P = polyfit(x_data,Y_data,1);
        y1_all(n_cell) = P(1);
        y0_all(n_cell) = P(2);
        fitmeans(n_cell) = P(2) + P(1)*num_trials/2;
    end
    
    figure;
    s1 = subplot(2,1,1);
    plot(fitmeans, '.'); axis padded
    ylabel('Mean response');
    s2 = subplot(2,1,2);
    plot(y1_all, '.'); axis padded;
    ylabel('Linear fit slope');
    xlabel('Sorted cells');
    linkaxes([s1, s2], 'x');
    
    figure; hold on;
    plt1 = plot(y1_all, fitmeans, '.');
    xlabel('Linear fit slope');
    ylabel('Mean cell response');
    f = gcf;
    line([0 0], f.Children.YLim, 'color', 'k')
    line(f.Children.XLim, [0 0], 'color', 'k')
    legend('Cells')
    uistack(plt1,'top');
    axis tight;
    title('Mean cell resp vs resp change with trials');
    
    plot_cells = [1, 2, num_cells-1, num_cells];
    
    for n_cell_idx = 1:numel(plot_cells)
        n_cell = plot_cells(n_cell_idx);
        Y_data = mean_resp_by_trial(n_cell,:);
        Y_fit = x_data*y1_all(n_cell) + y0_all(n_cell);
        figure; hold on;
        plot(x_data, Y_data, '.');
        plot(x_data, Y_fit);
        title(sprintf('Cell %d change in resp magnitude with trials', n_cell)); ylim(ylims1);
        xlabel('Trials');
        ylabel('Response magnitude')
        f = gcf;
        text(f.Children.XLim(2)-diff(f.Children.XLim)*.4, f.Children.YLim(2)-diff(f.Children.YLim)*.03, sprintf('Linear fit slope = %.2g', y1_all(n_cell)))
    end
 
    A_all = cell(ddata.num_planes,1);
    for n_pl = 1:ddata.num_planes
        A_all{n_pl} = full(ddata.OA_data{n_pl}.est.A);
    end

    color1 = jet(1000);

    A = cat(2,A_all{:});

    mean_mean_resp = mean(mean_resp);
    delta_resp = max(abs(mean_resp - mean_mean_resp));

    color_lut = linspace(mean_mean_resp - delta_resp, mean_mean_resp + delta_resp, 1000);

    [~, idx2] = min((mean_resp - color_lut).^2, [], 2);

    A_all2 = cell(ddata.num_planes,1);
    n_cell_all = 1;
    for n_pl = 1:ddata.num_planes
        A_flat = zeros(256 * 256, 3);
        for n_cell = 1:ddata.num_cells_pl{n_pl}
            temp_mask = logical(A_all{n_pl}(:,n_cell));
            A_flat(temp_mask,:) = A_all{n_pl}(temp_mask,n_cell).*color1(idx2(n_cell_all),:)./max(A_all{n_pl}(temp_mask,n_cell));
            n_cell_all = n_cell_all + 1;
        end
        A_flat2 = permute(reshape(A_flat, 256, 256, 3), [2, 1, 3]);
        if n_pl < ddata.num_planes
            A_all2{n_pl} = [A_flat2, ones(256, 1, 3)];
        else
            A_all2{n_pl} = A_flat2;
        end
    end

    A_all2 = cat(2, A_all2{:});

    figure; imagesc(A_all2); axis equal tight;
    title('opto response magnitudes');
    f1 = figure; 
    imagesc([], color_lut, reshape(color1, [], 1, 3))
    f1.Children.YDir = 'normal';
    ylabel('resp magnitude')
else
    disp('pockel chan has no stim')
end
end