function f_dv_dim_est_trial_shuff(app)

mean_win = [.1 .9];
samp_start = 5;
samp_interval = 5;
samp_max = 40;
dist_metric = app.DistmethodDropDown.Value; % pca isomap


[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = size(data,1);

tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';
num_tn = numel(tn_all);

params = f_dv_gather_params(app);


data_all = cell(num_dsets,3);
added_data = false(num_dsets,1);
for n_dset = 1:num_dsets
    params.n_dset = n_dset;
    ddata = data(n_dset,:);
    stats1 = cat(1,ddata.stats{n_pl});
    cdata = f_dv_compute_cdata(ddata, params);
   
    firing_rate = cat(1,cdata.S_sm);
    
    [sel_cells] = f_dv_get_resp_vals_cells(app, stats1, tn_all);
    sel_cells2 = logical(sum(sel_cells,2));
    
    num_cells = sum(sel_cells2);
    
    if num_cells>=10
        added_data(n_dset) = 1;
        
        firing_rate2 = firing_rate(sel_cells2,:);



        if isempty(tn_all)
            disp('Specified trial type does not exist');
        else
            stim_times = ddata.stim_frame_index{1};
            mmn_freq = ddata.MMN_freq{1};

            [~, trial_frames] = f_dv_compute_window_t(mean_win, cdata(1).volume_period);

            trial_types = ddata.trial_types{1};
            trial_data_sort = f_get_stim_trig_resp(firing_rate2, stim_times, trial_frames);

            if ~isempty(mmn_freq)
                [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);
            else
                trial_data_sort_wctx = trial_data_sort;
                trial_types_wctx = trial_types;
            end

            trial_data_sort_wctx_2d = reshape(mean(trial_data_sort_wctx,2), num_cells, []);

            trial_idx = logical(sum(tt_all == trial_types_wctx,2));
            trial_seq = trial_types_wctx(trial_idx);
            trial_data_sort2 = trial_data_sort_wctx_2d(:,trial_idx);

            tn_seq = f_tt_to_tn(trial_seq, app.ops,0);

    %         if app.sortbytrialtypeCheckBox.Value
    %             [tn_seq, sort_idx] = sort(tn_seq);
    %             trial_data_sort2 = trial_data_sort2(:,:,sort_idx);
    %             trial_seq = trial_seq(sort_idx);
    %         end
        end
   
        samps_num = samp_start:samp_interval:min(samp_max, num_cells);
        num_samp = sum(samps_num);
        num_samp_gr = numel(samps_num);

        dist_mean_all = zeros(num_tn, num_samp_gr);
        dist_std_all = zeros(num_tn, num_samp_gr);
        for n_tn = 1:numel(tn_all)
            idx1 = tn_seq == tn_all(n_tn);

            samps = randsample(num_cells, num_samp, 1);
            n_samp_start = 1;
            for n_samp_gr = 1:num_samp_gr

                samp1 = samps(n_samp_start:(n_samp_start+samps_num(n_samp_gr)-1));
                trial_data_sort3 = trial_data_sort2(samp1,idx1);

                vec_mean = mean(trial_data_sort3,2);

                dist1 = pdist2(vec_mean', trial_data_sort3', dist_metric);  % euclidean cosine

                dist_mean_all(n_tn, n_samp_gr) = mean(dist1);
                dist_std_all(n_tn, n_samp_gr) = std(dist1);

                n_samp_start = n_samp_start + samps_num(n_samp_gr)-1;
            end

        end

        data_all{n_dset,1} = samps_num;
        data_all{n_dset,2} = dist_mean_all;
        data_all{n_dset,3} = dist_std_all;
    end
end

data_all2 = data_all(added_data,:);

num_dsets2 = sum(added_data);

samps_num2 = samp_start:samp_interval:samp_max;
num_samp_gr = numel(samps_num2);

figure; hold on;
for n_tn = 1:num_tn
    color1 = app.ops.context_types_all_colors2{tn_all(n_tn)};
    dist_all4 = nan(num_dsets2, num_samp_gr);
    for n_dset = 1:num_dsets2
        data3 = data_all2{n_dset, 2}(n_tn,:);
        if ~sum(isnan(data3))
            dist_all4(n_dset,1:numel(data3)) = data3;
            plot(data_all2{n_dset, 1}, data3, 'color', [color1, 0.4])
        end
    end
    data_means = nan(1, num_samp_gr);
    for n_gr = 1:num_samp_gr
        data_means(n_gr) = mean(dist_all4(~isnan(dist_all4(:,n_gr)),n_gr));
    end
    idx1 = ~isnan(data_means);
    plot(samps_num2(idx1), data_means(idx1), 'Linewidth', 2, 'color', color1);
end
title('mean dist')

% figure; hold on;
% for n_tn = 1:num_tn
%     color1 = app.ops.context_types_all_colors2{tn_all(n_tn)};
%     dist_all4 = nan(num_dsets2, num_samp_gr);
%     for n_dset = 1:num_dsets2
%         data3 = data_all2{n_dset, 3}(n_tn,:);
%         if ~sum(isnan(data3))
%             dist_all4(n_dset,1:numel(data3)) = data3;
%             plot(data_all2{n_dset, 1}, data3, 'color', [color1, 0.4])
%         end
%     end
%     data_means = nan(1, num_samp_gr);
%     for n_gr = 1:num_samp_gr
%         data_means(n_gr) = mean(dist_all4(~isnan(dist_all4(:,n_gr)),n_gr));
%     end
%     idx1 = ~isnan(data_means);
%     plot(samps_num2(idx1), data_means(idx1), 'Linewidth', 2, 'color', color1);
% end
% title('dist std')


end