function f_dv_similarity_onevone(data0, params, ops)

%tn_all = f_dv_get_trial_number(params);
[data, title_tag] = f_dv_get_data_by_mouse_selection(data0, params);
%[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);

dist_metric = params.distance_method; % pca isomap

do_fdr = 1;

markSize = 4;

%[region_num, reg_tag, leg_list, reg_col] = f_dv_get_region_sel_val(app);
region_num = [1, 2, 3, 4];

shuff_rep = 100;

if ~strcmpi(params.responsive_cells_select, 'All')
    params.responsive_cells_select = 'Resp marg';
end

num_dsets = size(data,1);
%tn_all = f_dv_get_trial_number(app);
tn_all = [18, 19, 20, 28, 29, 30];
[num_gr, num_tn] = size(tn_all);
num_effdsets = num_dsets*num_gr;

%num_cont = 8;
%tn_all = [1:num_cont 18 19 20 28 29 30]; %f_dv_get_trial_number(app);

if params.do_similarity
    method_tag1 = sprintf('%s similarity', dist_metric);
else
    method_tag1 = sprintf('%s distance', dist_metric);
end
params.region = 'all comb';
[~, ~, ~, features_proc] = f_dv_get_feature(params.plot_feature, tn_all, data, params, ops);
features_proc2 = reshape(features_proc, num_effdsets, num_tn);

features_proc3 = cell(num_dsets,1);
for n_dset = 1:num_dsets
    features_proc3{n_dset} = cat(2,features_proc2{n_dset,:});
end

features_proc4 = cat(1, features_proc3{:});

hasnan1 = logical(sum(isnan(features_proc4),2));
features_proc5 = features_proc4(~hasnan1,:);

dist2 = squareform(pdist(features_proc5', dist_metric));
% 
% data_all = cell(num_dsets,1);
% cell_dset_idx = cell(num_dsets,1);
% num_cells_all = zeros(num_dsets,1);
% for n_dset = 1:num_dsets
%     stats1 = cat(1,data(n_dset,:).stats{n_pl});
% 
%     [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], resp_cell_sel);
% 
%     resp_vals2 = cat(2,resp_vals{:});
%     num_cells_all(n_dset) = size(resp_vals2,1);
%     cell_dset_idx{n_dset} = ones(num_cells_all(n_dset),1)*n_dset;
% 
%     data_all{n_dset} = resp_vals2;
% end
% 
% 
% data_all2 = cat(1,data_all{:});
% hasnan1 = logical(sum(isnan(data_all2),2));
% data_all2 = data_all2(~hasnan1,:);
% 
% num_cells = size(data_all2,1);
% 
% data_all2_shuff = zeros(num_cells, num_tn);
% for n_cell = 1:num_cells
%     idx_sh = randsample(6,6, false);
%     data_all2_shuff(n_cell,:) = data_all2(n_cell,idx_sh);
% end
% 
% cell_dset_idx2 = cat(1, cell_dset_idx{:});
% cell_dset_idx2 = cell_dset_idx2(~hasnan1);
% 
% dist2 = squareform(pdist(data_all2', dist_metric));
% 
% %dist2_sh = squareform(pdist(data_all2_shuff', dist_metric));
% 
dist3 = if_get_mat_plot_data(dist2, params.mat_tri);

num_cells = size(features_proc5,1);

if strcmpi(dist_metric, 'euclidean')
    dist3 = dist3./sqrt(num_cells)*sqrt(50);
    ylab = sprintf('normalized %s', method_tag1);
else
    ylab = sprintf('%s', method_tag1);
end

if params.do_similarity
    dist4 = 1-dist3;
else
    dist4 = dist3;
end

if strcmpi(params.mat_tri, 'UTri')
    dist4(logical(tril(dist4))) = nan;
elseif strcmpi(params.mat_tri, 'LTri')
    dist4(logical(triu(dist4))) = nan;
end

figure;
imagesc(dist4);
axis square;
set(gca,'xtick',1:num_tn,'xticklabel',ops.context_types_labels(tn_all));
set(gca,'ytick',1:num_tn,'yticklabel',ops.context_types_labels(tn_all));
colorbar;
colormap(params.colormap)
title(sprintf('full %s; %s', ylab, title_tag), 'interpreter', 'none');
num_pairs = nchoosek(numel(tn_all),2);

idx_data = find(triu(squareform(pdist((1:numel(tn_all))', 'euclidean')))>0);
[row,col] = ind2sub([numel(tn_all),numel(tn_all)], idx_data);
%pairs = [row, col];

dist_all = nan(num_dsets, num_pairs);
sq_all = nan(num_tn, num_tn, num_dsets);
dist_all_sh = nan(num_dsets, num_pairs, shuff_rep);
sq_all_sh = nan(num_tn, num_tn, num_dsets, shuff_rep);

feature3_proc = cell(num_dsets,1);
features3_proc_sh = cell(num_dsets,1);
for n_dset = 1:num_dsets
    data_d = features_proc3{n_dset};
    
    hasnan1 = logical(sum(isnan(data_d),2));
    data_d2 = data_d(~hasnan1,:);

    feature3_proc{n_dset} = data_d2;
    
    num_cells = size(data_d2,1);

    if num_cells > 5
        
        data_d2_sh = zeros(num_cells, num_tn, shuff_rep);
        
        for n_shuff = 1:shuff_rep
            
            

            for n_cell = 1:num_cells
                idx_sh = randsample(6,6, false);
                data_d2_sh(n_cell,:,n_shuff) = data_d2(n_cell,idx_sh);
            end
            
            [dist3_sh, dist3_sh_th] = if_get_dist(data_d2_sh(:,:,n_shuff), dist_metric, 0, params.mat_tri);
            
            if strcmpi(dist_metric, 'euclidean')
                dist3_sh = dist3_sh./sqrt(num_cells)*sqrt(50);
            end
            
            dist_all_sh(n_dset,:, n_shuff) = dist3_sh(idx_data);
            sq_all_sh(:,:, n_dset, n_shuff) = dist3_sh_th;

        end
    
        features3_proc_sh{n_dset} = data_d2_sh;

        
        [dist3, dist3_th] = if_get_dist(data_d2, dist_metric, 0, params.mat_tri);
        
        if strcmpi(dist_metric, 'euclidean')
            dist3 = dist3./sqrt(num_cells)*sqrt(50);
        end
        
        dist_all(n_dset,:) = dist3(idx_data);
        sq_all(:,:, n_dset) = dist3_th;
    
        
    end
end

hasnan1 = logical(sum(isnan(dist_all),2));
dist_all2 = dist_all(~hasnan1,:);
dist_all2_sh = dist_all_sh(~hasnan1,:,:);
sq_all2 = sq_all(:,:, ~hasnan1);
sq_all2_sh = sq_all_sh(:,:, ~hasnan1,:);

num_dsets2 = sum(~hasnan1);

% normalize euqlidean
if strcmpi(dist_metric, 'euclidean')
    ylab = sprintf('normalized %s', method_tag1);
else
    ylab = sprintf('%s', method_tag1);
end

if params.do_similarity
    dist_all3 = 1-dist_all2;
    dist_all3_sh = 1 - dist_all2_sh;
    sq_all3 = 1 - sq_all2;
else
    dist_all3 = dist_all2;
    dist_all3_sh = dist_all2_sh;
    sq_all3 = sq_all2;
end

labels1 = ops.context_types_labels_trim(tn_all);
legend1 = cell(num_pairs,1);
for n_pair = 1:num_pairs
    legend1{n_pair} = sprintf('%s-%s', labels1{row(n_pair)}, labels1{col(n_pair)});
end

sq_mean = mean(sq_all3,3);
if strcmpi(params.mat_tri, 'UTri')
    sq_mean(logical(tril(sq_mean))) = nan;
elseif strcmpi(params.mat_tri, 'LTri')
    sq_mean(logical(triu(sq_mean))) = nan;
end

figure;
imagesc(sq_mean);
axis square;
set(gca,'xtick',1:num_tn,'xticklabel',ops.context_types_labels(tn_all));
set(gca,'ytick',1:num_tn,'yticklabel',ops.context_types_labels(tn_all));
colorbar;
colormap(params.colormap);
title(sprintf('mean %s; %s', method_tag1, title_tag), 'interpreter', 'none');

x_noise1 = (rand(size(dist_all3))-0.5)/5 + (1:num_pairs);
figure;
hold on;
plot(x_noise1', dist_all3', '.', color=[.7 .7 .7], markersize=markSize);
plot(mean(dist_all3,1, 'omitnan'), '_', color='black', linewidth=1, markersize=10);
errorbar((1:num_pairs), mean(dist_all3,1, 'omitnan'), std(dist_all3, [], 1, 'omitnan')/sqrt(num_dsets-1), '.', color='black', linewidth=1, markersize=10);
xlim([0, num_pairs+1]);
ylabel(ylab);
set(gca,'xtick',1:num_pairs,'xticklabel',legend1);
title(sprintf('%s; %s; %s', title_tag, method_tag1, params.plot_feature), 'interpreter', 'none');

% anova stats between groups... too many combinations
if 0 
    legend2 = repmat(legend1', num_dsets2, 1);

    [p_tun, tbl_tun, stats_tun]  = anova1(dist_all3(:), legend2(:), 'off');
    title_tag4 = sprintf('%s; %s', method_tag1, title_tag);
    f_dv_plot_anova1(p_tun, tbl_tun, stats_tun, title_tag4, legend1);

end

shuf_mean2 = mean(dist_all3_sh,3);
shuf_mean = mean(shuf_mean2,1, 'omitnan');
shuf_sem = std(shuf_mean2,1)/sqrt(num_dsets2-1);

dat_mean = mean(dist_all3,1, 'omitnan');
dat_sem = std(dist_all3, [], 1, 'omitnan')/sqrt(num_dsets2-1);

linew = 1.2;
marks = 15;
figure;
hold on;
x_noise1 = (rand(num_dsets2, num_pairs)-0.5)/5 + (1:num_pairs);
plot(x_noise1', shuf_mean2', '.', color=[.6 .6 .6], markersize=markSize);
plot(shuf_mean, '_', color='black', linewidth=linew, markersize=10);
pls = errorbar((1:num_pairs), shuf_mean, shuf_sem, '.', color='black', linewidth=linew, markersize=marks);
x_noise1 = (rand(size(dist_all3))-0.5)/5 + (1:num_pairs);
plot(x_noise1', dist_all3', '.', color=[.4 0.7 0.8], markersize=markSize);
plot(dat_mean, '_', color=[0 0.4470 0.7410], linewidth=linew, markersize=10);
pld = errorbar((1:num_pairs), dat_mean, dat_sem, '.', color=[0 0.4470 0.7410], linewidth=linew, markersize=marks);
xlim([0, num_pairs+1]);
ylabel(ylab);
set(gca,'xtick',1:num_pairs,'xticklabel',legend1);
title(sprintf('%s; %s; %s', title_tag, method_tag1, params.plot_feature), 'interpreter', 'none');
legend([pld(1), pls(1)], {'data', 'shuffle'})

if strcmpi(ylab, 'cosine distance')
    yyaxis left
    axl = gca;
    ylim1 = axl.YLim;
    yyaxis right
    
    ylabel('cosine similarity');
    axr = gca;
    axr.YLim =[1-ylim1(2),1-ylim1(1)];
    axr.YDir = 'reverse';
    linkaxes([axl axr]);
end

n_cat = 13;
fprintf('%s, %s cat mean=%.3f, sem=%.3f\n', ylab, legend1{n_cat}, dat_mean(n_cat), dat_sem(n_cat))
n_cat = 4;
fprintf('%s, %s cat mean=%.3f, sem=%.3f\n', ylab, legend1{n_cat}, dat_mean(n_cat), dat_sem(n_cat))

p_all = zeros(num_pairs,1);
t_all = zeros(num_pairs,1);

for  n_pair = 1:num_pairs
    data_temp = dist_all3(:,n_pair);
    data_temp_shuff = shuf_mean2(:,n_pair);
    
    %dat_mean = mean(data_temp);
    %dat_sem = std(data_temp);
    
    %shuf_mean = mean(data_temp_shuff);
    %shuf_sem = std(data_temp_shuff);
    
    [~, p_val, ~, stats] = ttest(data_temp, data_temp_shuff);
    p_all(n_pair) = p_val;
    t_all(n_pair) = stats.tstat;
    
end

if do_fdr
    p_vals_adj = f_FDR_correction(p_all);
    p_tag = 'p_adj';
    fdr_tag = 'FDR adj ';
else
    p_vals_adj = p_all;
    p_tag = 'p';
    fdr_tag = '';
end

fprintf('%s; %s\n%s paired ttest vs shuffled\n', method_tag1, title_tag, fdr_tag);

for  n_pair = 1:num_pairs
    p_val = p_vals_adj(n_pair);
    if p_val>0.01
        sig_tag = '*';
    elseif p_val>0.001
        sig_tag = '**';
    elseif p_val<=0.001
        sig_tag = '***';
    end
    
    if dat_mean(n_pair) > shuf_mean(n_pair)
        fprintf('%s > shuff t=%.3f; %s=%.2e%s\n', legend1{n_pair}, t_all(n_pair), p_tag, p_val, sig_tag)
    else
        fprintf('%s < shuff t=%.3f; %s=%.2e%s\n', legend1{n_pair}, -t_all(n_pair), p_tag, p_val, sig_tag)
    end
end

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
% [cdata, ~] = f_dv_get_new_cdata_stats(ddata, params);
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

function [dist3, dist3_th] = if_get_dist(data, dist_metric, do_similarity, mat_plot_type)

dist2 = squareform(pdist(data', dist_metric));
dist2_tr = if_get_mat_plot_data(dist2, mat_plot_type);

if do_similarity
    dist3 = 1-dist2;
    dist3_th = 1-dist2_tr;
else
    dist3 = dist2;
    dist3_th = dist2_tr;
end

end


function data_out = if_get_mat_plot_data(data_in, mat_plot_type)

if strcmpi(mat_plot_type, 'LTri')
    data_out = tril(data_in);
elseif strcmpi(mat_plot_type, 'UTri')
    data_out = triu(data_in);
elseif strcmpi(mat_plot_type, 'Full')
    data_out = data_in;
end

end
