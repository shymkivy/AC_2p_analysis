function f_dv_similarity_onevone(app)

dist_metric = app.LDdistmethodDropDown.Value; % pca isomap

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = size(data,1);

tn_all = f_dv_get_trial_number(app);

%num_cont = 8;
%tn_all = [1:num_cont 18 19 20 28 29 30]; %f_dv_get_trial_number(app);

title_tag1 = sprintf('%s; dist %s', title_tag, dist_metric);


data_all = cell(num_dsets,1);
cell_dset_idx = cell(num_dsets,1);
num_cells_all = zeros(num_dsets,1);
for n_dset = 1:num_dsets
    stats1 = cat(1,data(n_dset,:).stats{n_pl});

    if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
        resp_cell_sel = 'All';
    else
        resp_cell_sel = 'Resp marg';
    end

    [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], resp_cell_sel);

    resp_vals2 = cat(2,resp_vals{:});
    num_cells_all(n_dset) = size(resp_vals2,1);
    cell_dset_idx{n_dset} = ones(num_cells_all(n_dset),1)*n_dset;
   
    data_all{n_dset} = resp_vals2;
end


data_all2 = cat(1,data_all{:});
hasnan1 = logical(sum(isnan(data_all2),2));
data_all2 = data_all2(~hasnan1,:);

cell_dset_idx2 = cat(1, cell_dset_idx{:});
cell_dset_idx2 = cell_dset_idx2(~hasnan1);

dist2 = squareform(pdist(data_all2', dist_metric));

dist3 = tril(dist2);

figure;
imagesc(dist3);
axis square
set(gca,'xtick',1:4,'xticklabel',app.ops.context_types_labels(tn_all))
set(gca,'ytick',1:4,'yticklabel',app.ops.context_types_labels(tn_all))
title(sprintf('full dist; %s', title_tag1))
colorbar;
colormap('gray')

n_pairs = nchoosek(numel(tn_all),2);

idx_data = find(triu(squareform(pdist([1:numel(tn_all)]', 'euclidean')))>0);
[row,col] = ind2sub([numel(tn_all),numel(tn_all)], idx_data);
pairs = [row, col];


dist_all = zeros(num_dsets,n_pairs);
sq_all = zeros(numel(tn_all), numel(tn_all), num_dsets);
for n_dset = 1:num_dsets
    idx1 = cell_dset_idx2 == n_dset;
    data_all3 = data_all2(idx1,:);
    
    dist2 = squareform(pdist(data_all3', dist_metric));
    dist_all(n_dset,:) = dist2(idx_data);

    sq_all(:,:, n_dset) = tril(dist2);
end

% normalize euqlidean
if strcmpi(dist_metric, 'euclidean')
    dist_all = dist_all./sqrt(num_cells_all)*sqrt(mean(num_cells_all));
    ylab = 'normalized euclidean distance';
else
    ylab = sprintf('%s distance', dist_metric);
end


labels1 = app.ops.context_types_labels(tn_all);
legend1 = cell(numel(row),1);
for n_pair = 1:numel(row)
    legend1{n_pair} = sprintf('%s-%s', labels1{row(n_pair)}, labels1{col(n_pair)});
end

figure;
imagesc(mean(sq_all,3));
axis square
set(gca,'xtick',1:4,'xticklabel',app.ops.context_types_labels(tn_all))
set(gca,'ytick',1:4,'yticklabel',app.ops.context_types_labels(tn_all))
title(sprintf('mean dist; %s', title_tag1))
colorbar;
colormap('gray')

x_noise1 = (rand(size(dist_all))-0.5)/5 + [1:n_pairs];
figure;
hold on;
plot(x_noise1', dist_all', '.', color=[.7 .7 .7])
plot(mean(dist_all,1), '_', color='black', linewidth=2, markersize=15)
errorbar([1:n_pairs], mean(dist_all,1), std(dist_all, [], 1)/sqrt(num_dsets-1), '.', color='black', linewidth=1, markersize=10);
title(sprintf('%s', title_tag1));
ylabel(ylab);
set(gca,'xtick',1:n_pairs,'xticklabel',legend1);



% 
% cell_dset_idx2 = cat(1,cell_dset_idx{:});
% MMN_freq = data.MMN_freq;
% 
% tn_d = [18, 19, 20; 28, 29, 30];
% MMMN_id = [2 2 2; 1 1 1];
% MMMN_id_opp = [1 1 1; 2 2 2];
% 
% [num_mmn, num_tnd] = size(tn_d);
% dist_all = cell(num_dsets, num_mmn, num_tnd);
% 
% for n_dset = 1:num_dsets
%     title_tag2 = sprintf('%s; dset%d; %dcells', title_tag1, data.idx(n_dset), num_cells_all(n_dset));
%     resp_vals = data_all{n_dset};
% 
%     cont_vec = resp_vals(:,1:num_cont);
% 
%     for n_mmn = 1:num_mmn
%         for n_tnd = 1:num_tnd
%             tn_idx = tn_all == tn_d(n_mmn, n_tnd);
%             vec1 = resp_vals(:,tn_idx);
%             dist_all{n_dset, n_mmn, n_tnd} = pdist2(cont_vec', vec1', dist_metric);
%         end
%     end
% end
% 
% dist_all_same = cell(num_dsets, n_mmn, n_tnd);
% dist_all_opp = cell(num_dsets, n_mmn, n_tnd);
% 
% pad = 2;
% for n_dset = 1:num_dsets
% 
%     MMN_freq1 = MMN_freq{n_dset};
%     for n_mmn = 1:num_mmn
%         for n_tnd = 1:num_tnd
%             MMN_freq2 = MMN_freq1(MMMN_id(n_mmn, n_tnd));
%             MMN_freq2_opp = MMN_freq1(MMMN_id_opp(n_mmn, n_tnd));
%             dist_temp = dist_all{n_dset, n_mmn, n_tnd};
% 
%             dist_temp2 = ones(pad*2+1, 1);
%             for n_pad = -pad:pad
%                 if and(MMN_freq2+n_pad >= 1 , MMN_freq2+n_pad <= num_cont)
%                     dist_temp2(pad+1+n_pad) = dist_temp(MMN_freq2+n_pad);
%                 end
%             end
%             dist_all_same{n_dset, n_mmn, n_tnd} = dist_temp2;
% 
% 
%             dist_temp2 = ones(pad*2+1, 1);
%             for n_pad = -pad:pad
%                 if and(MMN_freq2_opp+n_pad >= 1 , MMN_freq2_opp+n_pad <= num_cont)
%                     dist_temp2(pad+1+n_pad) = dist_temp(MMN_freq2_opp+n_pad);
%                 end
%             end
%             dist_all_opp{n_dset, n_mmn, n_tnd} = dist_temp2;
%         end
%     end
% end
% 
% 
% 
% x_lab = -pad:pad;
% figure; hold on;
% for n_tnd = 1:num_tnd
%     dist_temp = cat(2, dist_all_same{:, :, n_tnd});
%     dist_temp_opp = cat(2, dist_all_opp{:, :, n_tnd});
% 
%     plot(x_lab, dist_temp, '-', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)} 0.2])
%     plot(x_lab, mean(dist_temp,2), 'o-', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
% 
%     plot(x_lab, dist_temp_opp, '--', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 0.2])
%     plot(x_lab, mean(dist_temp_opp,2), 'o--', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
% end
% xlabel('Frequency around MMN');
% ylabel([dist_metric ' distance']);
% title(title_tag1)
% 
% 
% %%
% MMN_freq5 = cat(1, MMN_freq{:});
% 
% % first location in mmn is 2, second is 1
% high_idx = 2 - double(MMN_freq5(:,1) < MMN_freq5(:,2));
% low_idx = 3 - high_idx;
% 
% dist_temp_high_same = cell(num_dsets, n_tnd);
% dist_temp_high_opp = cell(num_dsets, n_tnd);
% dist_temp_low_same = cell(num_dsets, n_tnd);
% dist_temp_low_opp = cell(num_dsets, n_tnd);
% 
% for n_dset = 1:num_dsets
%     for n_tnd = 1:num_tnd
%         dist_temp_high_same{n_dset, n_tnd} = dist_all_same{n_dset, high_idx(n_dset), n_tnd};
%         dist_temp_high_opp{n_dset, n_tnd} = dist_all_opp{n_dset, high_idx(n_dset), n_tnd};
% 
%         dist_temp_low_same{n_dset, n_tnd} = dist_all_same{n_dset, low_idx(n_dset), n_tnd};
%         dist_temp_low_opp{n_dset, n_tnd} = dist_all_opp{n_dset, low_idx(n_dset), n_tnd};
%     end
% end
% 
% 
% x_lab = -pad:pad;
% figure; hold on;
% for n_tnd = 1:num_tnd
% 
%     dist_temp = cat(2, dist_temp_high_same{:, n_tnd});
%     dist_temp_opp = cat(2, dist_temp_high_opp{:, n_tnd});
% 
%     plot(x_lab, dist_temp, '-', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)} 0.2])
%     plot(x_lab, mean(dist_temp,2), 'o-', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
% 
%     plot(x_lab, dist_temp_opp, '--', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 0.2])
%     plot(x_lab, mean(dist_temp_opp,2), 'o--', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
% end
% xlabel('High frequency MMN');
% ylabel([dist_metric ' distance']);
% title([title_tag1 '; high MMN'])
% 
% x_lab = -pad:pad;
% figure; hold on;
% for n_tnd = 1:num_tnd
% 
%     dist_temp = cat(2, dist_temp_low_same{:, n_tnd});
%     dist_temp_opp = cat(2, dist_temp_low_opp{:, n_tnd});
% 
%     plot(x_lab, dist_temp, '-', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)} 0.2])
%     plot(x_lab, mean(dist_temp,2), 'o-', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
% 
%     plot(x_lab, dist_temp_opp, '--', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 0.2])
%     plot(x_lab, mean(dist_temp_opp,2), 'o--', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
% end
% xlabel('Low frequency MMN');
% ylabel([dist_metric ' distance']);
% title([title_tag1 '; low MMN'])

% 
% for n_dset = 1:num_dsets
%     figure; hold on;
%     for n_tnd = 1:num_tnd
%         if tn_d(n_tnd) > 20
%             pl_line = 'o--';
%         else
%             pl_line = 'o-';
%         end
%         plot(dist_all{n_dset, n_tnd}, pl_line, 'color', app.ops.context_types_all_colors2{tn_d(n_tnd)})
%     end
%     title(title_tag2, 'interpreter', 'none');
% end



%%

% decoder_type = 'bayes'; % tree, svm, bayes
% by_frame = 1;
% 
% tn_all = f_dv_get_trial_number(app);
% tt_all = app.ops.context_types_all(tn_all)';
% [num_gr, num_tn] = size(tn_all);
% 
% [data, title_tag] = f_dv_get_data_by_mouse_selection(app);
% num_dsets = size(data,1);
% params = f_dv_gather_params(app);
% 
% [region_num, reg_tag] = f_dv_get_region_sel_val(app);
% num_regions = numel(region_num);
% reg_all = app.ops.regions_to_analyze;
% 
% ddata = data(1,:);
% [cdata, ~] = f_dv_get_new_cdata_stats(app, ddata, params);
% 
% if by_frame
%     trial_window = [-1, 3];
%     [plot_t, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1).volume_period);
%     dec_acc_frames = zeros(num_dsets,sum(trial_frames), num_gr, num_regions);
%     dec_acc_frames_shuff = zeros(num_dsets,sum(trial_frames), num_gr, num_regions);
%     dec_acc_frames_bycl = zeros(num_dsets,sum(trial_frames), num_gr, num_regions, num_tn);
%     dec_acc_frames_bycl_shuff = zeros(num_dsets,sum(trial_frames), num_gr, num_regions, num_tn);
% else
%     trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
%     [plot_t, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1).volume_period);
%     dec_acc_curr_last_shuff = zeros(num_dsets,3, num_gr, num_regions);
% end
% 
% title_tag2 = sprintf('%s; %s reg; %dms smooth', title_tag, reg_tag, app.SmoothsigmamsEditField.Value);
% 
% mouse_id = cell(num_gr, num_dsets);
% dset_id = zeros(num_gr, num_dsets);
% num_cells_avail = zeros(num_gr, num_dsets, num_regions);
% done_dec = false(num_gr, num_dsets, num_regions);
% 
% fprintf('dset n/%d: ', num_dsets);
% 
% for n_dset = 1:num_dsets
% 
%     fprintf('%d..', n_dset);
% 
%     ddata = data(n_dset,:);
%     [cdata, stats1] = f_dv_get_new_cdata_stats(app, ddata, params);
%     mmn_freq = ddata.MMN_freq{1};
%     stim_times = ddata.stim_frame_index{1};
% 
%     firing_rate = cat(1,cdata.S_sm);
% 
%     trial_types = ddata.trial_types{1};
% 
% 
%     num_cells = sum([stats1.num_cells]);
% 
%     if app.UseregdatalabelsCheckBox.Value
%         if ~isempty(ddata.registered_data{1})
%             reg_cell_labels = ddata.registered_data{1}.reg_labels;
%         else
%             reg_cell_labels = zeros(num_cells,1);
%         end
%     else
%         reg_idx = find(strcmpi(reg_all, ddata.area));
%         reg_cell_labels = ones(num_cells,1)*reg_idx;
%     end
% 
%     for n_gr = 1:num_gr
%         tn1 = tn_all(n_gr, :);
%         mouse_id{n_gr, n_dset} = ddata.mouse_id{1};
%         dset_id(n_gr, n_dset) = n_dset;
% 
%         resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn1);
%         resp_cells2 = logical(sum(resp_cells,2));
% 
% 
%         for n_reg = 1:num_regions
%             reg_idx = reg_cell_labels == region_num(n_reg);
%             resp_reg_cell = and(resp_cells2, reg_idx);
% 
%             num_cells_avail(n_gr, n_dset, n_reg) = sum(resp_reg_cell);
% 
% 
%             stats1.peak_in_resp_win
% 
%             if num_cells_avail(n_gr, n_dset, n_reg) > 10
% 
%                 done_dec(n_gr, n_dset, n_reg) = 1;
% 
%                 firing_rate2 = firing_rate(resp_reg_cell,:);
% 
%                 trial_data_sort = f_get_stim_trig_resp(firing_rate2, stim_times, trial_frames);
% 
%                 trial_types2 = trial_types(2:400);
%                 trial_types2_last = trial_types(1:399);
%                 trial_types2_shuff = trial_types2(randperm(numel(trial_types2)));
% 
%                 trial_data_sort2 = trial_data_sort(:,:,2:400);
% 
%                 if by_frame
%                     for n_fr = 1:sum(sum(trial_frames))
% 
%                         trial_data_sort3 = squeeze(trial_data_sort2(:,n_fr,:));
% 
%                         trial_data_sort3n = trial_data_sort3 - mean(trial_data_sort3(:));
%                         trial_data_sort3n = trial_data_sort3n/std(trial_data_sort3n(:));
% 
%                         if strcmpi(decoder_type, 'tree')
%                             [trainedClassifier_tree, dec_acc_frames(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl(n_dset,n_fr,n_gr,n_reg,:)] = decoder_tree(trial_data_sort3n', trial_types2);
%                             [trainedClassifier_tree_shuff, dec_acc_frames_shuff(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl_shuff(n_dset,n_fr,n_gr,n_reg,:)] = decoder_tree(trial_data_sort3n', trial_types2_shuff);
%                         elseif strcmpi(decoder_type, 'bayes')
%                             [trainedClassifier_bayes, dec_acc_frames(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl(n_dset,n_fr,n_gr,n_reg,:)] = decoder_naivebayes(trial_data_sort3n', trial_types2);
%                             [trainedClassifier_bayes_shuff, dec_acc_frames_shuff(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl_shuff(n_dset,n_fr,n_gr,n_reg,:)] = decoder_naivebayes(trial_data_sort3n', trial_types2_shuff);
%                         elseif strcmpi(decoder_type, 'svm')
%                             [trainedClassifier_svm, dec_acc_frames(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl(n_dset,n_fr,n_gr,n_reg,:)] = decoder_svm(trial_data_sort3n', trial_types2, 1);
%                             [trainedClassifier_svm_shuff, dec_acc_frames_shuff(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl_shuff(n_dset,n_fr,n_gr,n_reg,:)] = decoder_svm(trial_data_sort3n', trial_types2_shuff, 1);
%                         end
%                     end
%                 else
%                     trial_data_sort3 = squeeze(mean(trial_data_sort2,2));
% 
%                     trial_data_sort3n = trial_data_sort3 - mean(trial_data_sort3(:));
%                     trial_data_sort3n = trial_data_sort3n/std(trial_data_sort3n(:));
% 
%                     if sprintf(decoder_type, 'tree')
%                         [trainedClassifier_tree, dec_acc_curr_last_shuff(n_dset,1,n_gr,n_reg)] = decoder_tree(trial_data_sort3n', trial_types2);
%                         [trainedClassifier_tree_last, dec_acc_curr_last_shuff(n_dset,2,n_gr,n_reg)] = decoder_tree(trial_data_sort3n', trial_types2_last);
%                         [trainedClassifier_tree_shuff, dec_acc_curr_last_shuff(n_dset,3,n_gr,n_reg)] = decoder_tree(trial_data_sort3n', trial_types2_shuff);
%                     elseif sprintf(decoder_type, 'bayes')
%                         [trainedClassifier_bayes, dec_acc_curr_last_shuff(n_dset,1,n_gr,n_reg)] = decoder_naivebayes(trial_data_sort3n', trial_types2);
%                         [trainedClassifier_bayes_last, dec_acc_curr_last_shuff(n_dset,2,n_gr,n_reg)] = decoder_naivebayes(trial_data_sort3n', trial_types2_last);
%                         [trainedClassifier_bayes_shuff, dec_acc_curr_last_shuff(n_dset,3,n_gr,n_reg)] = decoder_naivebayes(trial_data_sort3n', trial_types2_shuff);
%                     elseif sprintf(decoder_type, 'svm')
%                         [trainedClassifier_svm, dec_acc_curr_last_shuff(n_dset,1,n_gr,n_reg)] = decoder_svm(trial_data_sort3n', trial_types2, 1);
%                         [trainedClassifier_svm_last, dec_acc_curr_last_shuff(n_dset,2,n_gr,n_reg)] = decoder_svm(trial_data_sort3n', trial_types2_last, 1);
%                         [trainedClassifier_svm_shuff, dec_acc_curr_last_shuff(n_dset,3,n_gr,n_reg)] = decoder_svm(trial_data_sort3n', trial_types2_shuff, 1);
%                     end
%                 end
%                 % % mv regression
%                 % MnrModel_cont = fitmnr(trial_data_sort3n', trial_types2);
%                 % 
%                 % 
%                 % [d,p,stats] = manova1(trial_data_sort3n', trial_types2);
%                 % 
%                 % load fisheriris
%                 % 
%                 % MnrModel = fitmnr(meas,species);
%                 % MnrModel.Coefficients
% 
%                 % beta = mvregress(X,Y)
%             end
%         end
%     end  
% end
% 
% fprintf('\n');
% 
% 
% 
% for n_gr = 1:num_gr
%     for n_reg = 1:num_regions
%         done_dec2 = done_dec(n_gr, :, n_reg);
%         dec_acc_frames2 = dec_acc_frames(done_dec2,:,n_gr,n_reg);
%         dec_acc_frames_shuff2 = dec_acc_frames_shuff(done_dec2,:,n_gr,n_reg);
%         if by_frame
%             if ~isempty(dec_acc_frames_shuff2)
%                 figure; hold on; axis tight
%                 plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
%                 plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2)
% 
%                 plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
%                 plot(plot_t, mean(dec_acc_frames2,1), color=[0 0.4470 0.7410], LineWidth=2)
%                 title(sprintf('%s; %s decoder, freqs; %s', reg_all{region_num(n_reg)}, decoder_type, title_tag2), 'interpreter', 'none')
%             end
%         else
%             mean(dec_acc_curr_last_shuff)
%         end
%     end
% end
% 
% for n_gr = 1:num_gr
%     figure; hold on; axis tight
%     pl_all = cell(num_regions+1,1);
%     has_reg_data = false(num_regions+1,1);
%     for n_reg = 1:num_regions
%         done_dec2 = done_dec(n_gr, :, n_reg);
%         dec_acc_frames2 = dec_acc_frames(done_dec2,:,n_gr,n_reg);
%         dec_acc_frames_shuff2 = dec_acc_frames_shuff(done_dec2,:,n_gr,n_reg);
%         if by_frame
%             if ~isempty(dec_acc_frames_shuff2)
%                 has_reg_data(n_reg) = 1;
%                 has_reg_data(num_regions+1) = 1;
%                 %plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
%                 pl_all{num_regions+1} = plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2);
%                 col2 = app.ops.cond_colors{region_num(n_reg)};
%                 %plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
%                 pl_all{n_reg} = plot(plot_t, mean(dec_acc_frames2,1), color=col2, LineWidth=2);
%             end
%         else
%             mean(dec_acc_curr_last_shuff)
%         end
%     end
%     title(sprintf('%s decoder, freqs; %s', decoder_type, title_tag2), 'interpreter', 'none')
%     legend([pl_all{has_reg_data}], [reg_all(region_num(has_reg_data(1:num_regions))); {'Shuffle'}])
% end
% 
% for n_gr = 1:num_gr
%     for n_reg = 1:num_regions
%         done_dec2 = done_dec(n_gr, :, n_reg);
%         dec_acc_frames_bycl2 = dec_acc_frames_bycl(done_dec2,:,n_gr,n_reg,:);
%         dec_acc_frames_bycl_shuff2 = dec_acc_frames_bycl_shuff(done_dec2,:,n_gr,n_reg,:);
%         if by_frame
%             if ~isempty(dec_acc_frames_bycl2)
%                 figure; hold on; axis tight
%                 pl_all = cell(num_tn+1);
%                 for n_tt = 1:num_tn
%                     %plot(plot_t, dec_acc_frames_bycl_shuff2', color=[0 0 0 0.2])
%                     pl_all{num_tn+1} = plot(plot_t, mean(dec_acc_frames_bycl_shuff2(:,:,:,:,n_tt),1), color=[0 0 0], LineWidth=2);
%                 end
%                 for n_tt = 1:num_tn
%                     col2 = app.ops.context_types_all_colors2{tn_all(n_tt)};
%                     %plot(plot_t, dec_acc_frames_bycl2', color=[0 0.4470 0.7410 0.2])
%                     pl_all{n_tt} = plot(plot_t, mean(dec_acc_frames_bycl2(:,:,:,:,n_tt),1), color=col2, LineWidth=2);
%                 end
%                 title(sprintf('%s; %s decoder, freqs; %s', reg_all{region_num(n_reg)}, decoder_type, title_tag2), 'interpreter', 'none');
%                 legend([pl_all{:}], [app.ops.context_types_labels(tn_all); {'Shuffle'}]);
%             end
%         else
%             mean(dec_acc_curr_last_shuff)
%         end
%     end
% end

end