function f_dv_plot_raster(app)
%% plot raster of firing rates without enseble analysis

tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';

ddata = app.ddata;
cdata = f_dv_get_cdata(app);
firing_rate = cat(1,cdata.S_sm);

stats1 = cat(1,ddata.stats{:});

resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all);
resp_cells2 = logical(sum(resp_cells,2));

firing_rate2 = firing_rate(resp_cells2,:);

num_cells = sum(resp_cells2);

if app.shufflecellsCheckBox.Value
    firing_rate2 = firing_rate2(randperm(num_cells),:);
end

if isempty(tn_all)
    disp('Analyzing full trace')
    firing_rate3 = firing_rate2;
    tn_seq_plot = [];
else
    if isempty(tt_all)
        disp('Specified trial type does not exist');
    else
        stim_times = ddata.stim_frame_index{1};
        mmn_freq = ddata.MMN_freq{1};
        
        trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
        [~, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1).volume_period);

        trial_types = app.ddata.trial_types{1};
        trial_data_sort = f_get_stim_trig_resp(firing_rate2, stim_times, trial_frames);
        
        num_t = size(trial_data_sort,2);
        
        trial_types_ctx2 = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);
        trial_types_all = [trial_types, trial_types_ctx2];
        tr_last = [0; trial_types(1:end-1)];% .* double(trial_types<=10);

        idx2 = logical(sum(reshape(tt_all, [1, 1, numel(tt_all)]) == trial_types_all,3));

        trial_seq = cat(1, trial_types(idx2(:,1)), trial_types_ctx2(idx2(:,2),1), trial_types_ctx2(idx2(:,3),2));
        trial_data_sort2 = cat(3, trial_data_sort(:,:,idx2(:,1)), trial_data_sort(:,:,idx2(:,2)), trial_data_sort(:,:,idx2(:,3)));
        tr_last2 = cat(1, tr_last(idx2(:,1)), tr_last(idx2(:,2)), tr_last(idx2(:,3)));
        

        % if ~isempty(mmn_freq)
        %     [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);
        % else
        %     trial_data_sort_wctx = trial_data_sort;
        %     trial_types_wctx = trial_types;
        % end
        % trial_idx = logical(sum(tt_all == trial_types_wctx,2));
        % trial_seq = trial_types_wctx(trial_idx);
        % trial_data_sort2 = trial_data_sort_wctx(:,:,trial_idx);

        if app.shuffletrialsCheckBox.Value
            num_trials = numel(trial_seq);
            shuff_seq = randperm(num_trials);
            trial_data_sort2 = trial_data_sort2(:,:,shuff_seq);
            trial_seq = trial_seq(shuff_seq);
            tr_last2 = tr_last2(shuff_seq);
        end

        tn_seq = f_tt_to_tn(trial_seq, app.ops,0);
        tn_last_seq = f_tt_to_tn(tr_last2, app.ops, 1);

        if app.sortbytrialtypeCheckBox.Value
            if app.sortprevtrialCheckBox.Value
                [tn_last_seq, tr_last_sort_idx] = sort(tn_last_seq);
                tn_seq = tn_seq(tr_last_sort_idx);
                trial_data_sort2 = trial_data_sort2(:,:,tr_last_sort_idx);
                trial_seq = trial_seq(tr_last_sort_idx);
            end

            [tn_seq, sort_idx] = sort(tn_seq);
            trial_data_sort2 = trial_data_sort2(:,:,sort_idx);
            trial_seq = trial_seq(sort_idx);
            tn_last_seq = tn_last_seq(sort_idx);
            
        end
        
        tn_seq_plot = reshape(repmat(tn_seq, [1 num_t])',[],1);
        firing_rate3 = reshape(trial_data_sort2, num_cells, []);
        if app.plotprevtrialCheckBox.Value
            tn_seq_plot = [reshape(repmat(tn_last_seq, [1 num_t])',[],1), tn_seq_plot];
        end

        %% clustering and plotting
        if app.sortcellbytraceaveCheckBox.Value
            sort_data = squeeze(mean(trial_data_sort2,2));
        else
            sort_data = firing_rate3;
        end

        % remove inactive cells
        active_cells = sum(sort_data,2) ~= 0;
        sort_data(~active_cells,:) = [];
        firing_rate3(~active_cells,:) = [];
        
        num_cells = sum(active_cells);

        if app.sortbycellsimilarityCheckBox.Value
            hc_params.method = 'ward'; % ward(inner square), average, single(shortest)
            hc_params.distance_metric = 'cosine'; % none, euclidean, squaredeuclidean, cosine, hammilarity
            hc_params.plot_dist_mat = app.plotsortingstuffCheckBox.Value;
            hc_params.plot_clusters = app.plotsortingstuffCheckBox.Value;
            hc_params.num_clust = [];
            hc_params.title_tag = 'Coeffs (cells)';
            hclust_out_cell = f_hcluster_wrap(sort_data, hc_params);
            ord_cell = hclust_out_cell.dend_order;
        else
            ord_cell = 1:num_cells;
        end
        
        if ~isempty(mmn_freq)
            title_tag2 = sprintf('%s raster; mmn %d %d', app.ddata.dset_name_full{1}, mmn_freq(1), mmn_freq(2));
        else
            title_tag2 = sprintf('%s raster', app.ddata.dset_name_full{1});
        end
        
        vol_per = mean(cat(1, cdata.volume_period));

        num_frames = size(firing_rate3,2);
        raster_t = (1:num_frames)/(1000/vol_per);
        
        f_plot_raster_mean(firing_rate3(ord_cell,:), 1, raster_t, tn_seq_plot, app.ops.context_types_all_colors2, app.ColormapDropDown.Value, app.InvertcmapCheckBox.Value);
        sgtitle(title_tag2, 'interpreter', 'none');
        % clim1 = clim;
        % clim([0 clim1(2)]);

        if app.plotpopvecCheckBox.Value
            
            trial_data_sort2_mean = mean(trial_data_sort2,2);
            num_tn = numel(tt_all);
            trial_ave1 = zeros(num_cells, num_tn);
        
            for n_tn = 1:num_tn
                idx1 = trial_seq == tt_all(n_tn);
                trial_ave1(:,n_tn) = mean(trial_data_sort2_mean(:,:,idx1),3);
            end
            
            trial_ave1n = trial_ave1 - min(trial_ave1(:));
            trial_ave1n = trial_ave1n/max(trial_ave1n(:));
        
            tt_cont = tt_all(logical(sum(tt_all' == 1:10,2)));
            f1 = figure; hold on;
            for n_tn = 1:numel(tt_cont)
                idx2 = tt_all == tt_cont(n_tn);
                idx3 = app.ops.context_types_all == tt_cont(n_tn);
                plot(ones(num_cells,1)*n_tn + trial_ave1n(ord_cell, idx2), 1:num_cells, color=app.ops.context_types_all_colors2{idx3}, LineWidth=2)
                plot(ones(num_cells,1)*n_tn, 1:num_cells, color=[0.3 0.3 0.3 0.5], LineWidth=1, LineStyle=':')
            end
            f1.Children.YDir = 'reverse';
            axis tight
            xlim([.9, numel(tt_cont)+1]);
            title(sprintf('%s; cont vec', title_tag2), 'interpreter', 'none')

            tt_mmn = tt_all(logical(sum(tt_all' == app.ops.context_types_all([18, 19,20, 28, 29, 30])',2)));
            f1 = figure; hold on;
            for n_tn = 1:numel(tt_mmn)
                idx2 = tt_all == tt_mmn(n_tn);
                idx3 = find(app.ops.context_types_all == tt_mmn(n_tn),1);
                plot(ones(num_cells,1)*n_tn + trial_ave1n(ord_cell, idx2), 1:num_cells, color=app.ops.context_types_all_colors2{idx3}, LineWidth=2)
                plot(ones(num_cells,1)*n_tn, 1:num_cells, color=[0.3 0.3 0.3 0.5], LineWidth=1, LineStyle=':')
            end
            f1.Children.YDir = 'reverse';
            axis tight
            xlim([.9, numel(tt_mmn)+1]);
            title(sprintf('%s; ctx vec', title_tag2), 'interpreter', 'none')

        end
    end
end

end