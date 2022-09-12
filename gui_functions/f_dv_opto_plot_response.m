function f_dv_opto_plot_response(app)

n_pl = app.mplSpinner.Value;

ddata = app.ddata;

cdata = f_dv_get_cdata(app);
firing_rate = cat(1,cdata.S_sm);

base_resp_win = [5 20];
chan_idx = strcmpi(ddata.proc_ops{1}.chan_labels, 'Pockel');
stim_times = ddata.proc_data{1}.stim_times_frame{chan_idx, n_pl};

if ~isempty(stim_times)

    [~, T] = size(firing_rate);

    stim_times(stim_times<base_resp_win(1)) = [];
    stim_times(stim_times>(T-base_resp_win(2))) = [];

    sorted_data = f_get_stim_trig_resp(firing_rate, stim_times, base_resp_win);

    trial_ave = mean(sorted_data,3);
    mean_resp = mean((trial_ave(:,(base_resp_win(1)+1):(base_resp_win(1)+6))),2);

    [~, idx1] = sort(mean_resp, 'descend');

    figure; imagesc(trial_ave(idx1,:));
    xlabel('frames'); ylabel('sorted cells');
    title('opto triggered response')

    figure; hold on;
    plot(squeeze(sorted_data(idx1(2),:,:)), 'color', [.6 .6 .6]);
    plot(mean(sorted_data(idx1(1),:,:),3), 'm', 'linewidth', 2);
    axis tight; xlabel('frames'); ylabel('resp mag');
    title('top activated cell')

    figure; hold on;
    plot(squeeze(sorted_data(idx1(end),:,:)), 'color', [.6 .6 .6]);
    plot(mean(sorted_data(idx1(end),:,:),3), 'm', 'linewidth', 2);
    axis tight; xlabel('frames'); ylabel('resp mag');
    title('top inhibited cell')

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