function f_plot_stim_vec_trial_dist(trial_mean_vec, trial_raster, ops, dist_metric)

num_tt = size(trial_mean_vec{1},1);

dist_list = cell(numel(ops.regions_to_analyze),num_tt);
dist_mean_dset = cell(numel(ops.regions_to_analyze),num_tt);
for n_tt = 1:num_tt
    for n_cond = 1:numel(ops.regions_to_analyze)
        num_dsets = size(trial_mean_vec{n_cond},2);
        num_stim = numel(trial_raster{n_cond}{n_tt,1});
        
        dist_list2 = cell(num_stim, 1);
        dist_mean2 = zeros(num_stim, num_dsets);
        for n_stim = 1:num_stim
            dist_list3 = cell(num_dsets,1);
            for n_dset = 1:num_dsets
                temp_mean_vec = trial_mean_vec{n_cond}{n_tt,n_dset};
                temp_raster = trial_raster{n_cond}{n_tt, n_dset};
            
                dist_list3{n_dset} = pdist2(temp_mean_vec(:,n_stim)',temp_raster{n_stim}',dist_metric);
                dist_mean2(n_stim,n_dset) = mean(dist_list3{n_dset});
            end
            dist_list2{n_stim} = cat(2, dist_list3{:});
        end
        dist_list{n_cond, n_tt} = cat(2,dist_list2);
        dist_mean_dset{n_cond, n_tt} = dist_mean2;
        
        if strcmpi(ops.dred_params.trial_types_for_dist{n_tt}, 'mmn12')
            stim_ind = [1 4; 2 5; 3 6];
            figure; hold on;
            for n_stim = 1:size(stim_ind,1)
                [f, x] = ksdensity(cat(2,dist_list{n_cond, n_tt}{stim_ind(n_stim,1),:},dist_list{n_cond, n_tt}{stim_ind(n_stim,2),:}));
                plot(x, f/sum(f),'LineWidth', 2, 'color', ops.context_colors{n_stim});

                %h1 = histogram(cat(2,dist_list2{stim_ind(n_st,1),:},dist_list2{stim_ind(n_st,2),:}));
                %h1.Normalization = 'probability';
            end
            title([ops.regions_to_analyze{n_cond} ' trial ' dist_metric ' dist distribution']);
            ylim([0 .03]);
        end
    end
    
    if strcmpi(ops.dred_params.trial_types_for_dist{n_tt}, 'mmn12')
        stim_ind = [1 4; 2 5; 3 6];
        means1 = zeros(numel(ops.regions_to_analyze), size(stim_ind,1));
        sems1 = zeros(numel(ops.regions_to_analyze), size(stim_ind,1));
        for n_cond = 1:numel(ops.regions_to_analyze)
            
            for n_stim = 1:size(stim_ind,1)
                temp_list = cat(2,dist_list{n_cond, n_tt}{stim_ind(n_stim,1),:},dist_list{n_cond, n_tt}{stim_ind(n_stim,2),:});
                means1(n_cond,n_stim) = nanmean(temp_list);
                sems1(n_cond,n_stim) = nanstd(temp_list)/sqrt(numel(temp_list)-1);
            end
        end
        figure; hold on;
        for n_cond = 1:numel(ops.regions_to_analyze) 
            shadedErrorBar_YS(1:3, means1(n_cond,:), sems1(n_cond,:), ops.cond_colors{n_cond})
        end
        figure; hold on;
        b1 = bar(categorical({'Cont', 'Red', 'Dev'},{'Cont', 'Red', 'Dev'}),means1');
        for n_cond = 1:numel(ops.regions_to_analyze) 
            b1(n_cond).FaceColor = ops.cond_colors{n_cond};
        end
        ylim([0.2 0.6]);
        figure; hold on;
        curr_bar_ind = 1;
        for n_dd = 1:3
            for n_cond = 1:numel(ops.regions_to_analyze)
                b1 = bar(curr_bar_ind,means1(n_cond, n_dd));
                b1.FaceColor = ops.cond_colors{n_cond};
                er = errorbar(curr_bar_ind, means1(n_cond, n_dd), sems1(n_cond, n_dd));
                %er.Color = ops.cond_colors{n_cond};
                curr_bar_ind = curr_bar_ind + 1;
            end
            curr_bar_ind = curr_bar_ind + 1;
        end
        ylim([0.3 0.6]);
        title(['Cont Red Dev cosine distances'])
        
    else
        means1 = zeros(numel(ops.regions_to_analyze), numel(dist_list{n_cond, n_tt}));
        sems1 = zeros(numel(ops.regions_to_analyze), numel(dist_list{n_cond, n_tt}));
        for n_cond = 1:numel(ops.regions_to_analyze)       
            for n_stim = 1:numel(dist_list{n_cond, n_tt})
                means1(n_cond,n_stim) = nanmean(dist_list{n_cond, n_tt}{n_stim});
                sems1(n_cond,n_stim) = nanstd(dist_list{n_cond, n_tt}{n_stim})/sqrt(numel(dist_list{n_cond, n_tt}{n_stim})-1);
            end
        end
        figure; hold on;
        for n_cond = 1:numel(ops.regions_to_analyze) 
            shadedErrorBar_YS(1:10, means1(n_cond,:), sems1(n_cond,:), ops.cond_colors{n_cond})
        end
    end
%     figure; hold on;
%     b1 = bar(categorical(ops.regions_to_analyze,ops.regions_to_analyze), means1);
%     errorbar(categorical(ops.regions_to_analyze,ops.regions_to_analyze), means1,sems1)
%     for n_stim = 1:numel(dist_list{n_cond, n_tt})
%         b1(n_stim).FaceColor = ops.context_types_all_colors(n_stim,:);
%     end
%     ax1 = gca;
end

p_val_mat = f_get_tt_stats(corr_list);
figure;
subplot(2,1,2);
imagesc(p_val_mat);
caxis([0 1]);
sp1 = subplot(2,1,1);
f_plot_dset_deets(corr_list, ops, sp1)





end