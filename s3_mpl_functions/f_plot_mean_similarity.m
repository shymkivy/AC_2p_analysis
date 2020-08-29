function f_plot_mean_similarity(data, ops)

%% collect data
tt_types = fields(data);
tt_types_data = cell(numel(tt_types),1);
for n_tt = 1:numel(tt_types)
    tt_dist_cond = cell(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        chdata = data.(tt_types{n_tt}).(cond_name);
        
        num_dsets = numel(chdata.hclust_out_tr);
        dist1 = zeros(num_dsets,1);
        no_data = false(num_dsets,1);
        for n_dset = 1:num_dsets
            if ~isempty(chdata.hclust_out_tr{n_dset})
                dist1(n_dset) = nanmean(1-chdata.hclust_out_tr{n_dset}.dist);
            else
                no_data(n_dset) = true;
            end
        end
        tt_dist_cond{n_cond} = dist1(~no_data);
    end
    tt_types_data{n_tt} = tt_dist_cond;
end

%% create summed cond
combine_conds = {'dd1', 'dd2', 'dd1+2';...
                'cont1', 'cont2', 'cont1+2';...
                'red1', 'red2', 'red1+2'};

for n_comb = 1:size(combine_conds,1)
    if sum(strcmpi(tt_types, combine_conds{n_comb,1})) && sum(strcmpi(tt_types, combine_conds{n_comb,2}))
        size1 = numel(tt_types_data)+1;
        for n_cond = 1:numel(ops.regions_to_analyze)
            tt_types_data{size1,1}{n_cond,1} = [tt_types_data{strcmpi(tt_types, combine_conds{n_comb,1})}{n_cond}; tt_types_data{strcmpi(tt_types, combine_conds{n_comb,2})}{n_cond}];
        end 
        tt_types = cat(1, tt_types, combine_conds{n_comb,3});
    end
end

%% statistics?
% for n_tt =3:numel(tt_types)
%     mean_sims = cat(1,tt_types_data{n_tt}{2:3});
%     tt_types_list = cell(numel(tt_types_data{n_tt}),1);
%     for n_cond = 2:3%numel(tt_types_data{n_tt})
%         tt_types_list{n_cond} = repmat(ops.regions_to_analyze(n_cond),numel(tt_types_data{n_tt}{n_cond}),1);
%     end
%     mean_sims_labels = cat(1,tt_types_list{:});
%     
%     [~,~,stats]  = anova1(mean_sims, mean_sims_labels)
%     [c,~,~,gnames] = multcompare(stats)
%     [~, p] = ttest2(tt_types_data{n_tt}{2}, tt_types_data{n_tt}{3})
% end



%% plot
for n_tt = 1:numel(tt_types)
    if sum(tt_types{n_tt} == '+')
        figure;
        sp1 = subplot(2,1,1);
        f_plot_dset_deets(tt_types_data{n_tt}, ops, sp1);
        title(sprintf('trial-trial mean similarity, %s', tt_types{n_tt}));
        subplot(2,1,2)
        p_val_mat = f_get_tt_stats(tt_types_data{n_tt});
        imagesc(p_val_mat);
        caxis([0 1]);
    end
end

end