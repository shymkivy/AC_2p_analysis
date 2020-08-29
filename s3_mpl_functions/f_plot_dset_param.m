function f_plot_dset_param(data, param_name, ops)

%% collect data
tt_types = fields(data);
tt_types_data = cell(numel(tt_types),1);
for n_tt = 1:numel(tt_types)
    tt_cond_data = cell(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        chdata = data.(tt_types{n_tt}).(cond_name);
        num_dsets = numel(chdata.hclust_out_tr);
        tt_data = zeros(num_dsets,1);
        no_data = false(num_dsets,1);
        for n_dset = 1:num_dsets
            if ~isempty(chdata.hclust_out_tr{n_dset})
                tt_data(n_dset) = chdata.data_dim_est_full{n_dset}.(param_name);
            else
                no_data(n_dset) = true;
            end
        end
        tt_cond_data{n_cond} = tt_data(~no_data);
    end
    tt_types_data{n_tt,1} = tt_cond_data;
end

%% combine conds
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

%% plot
for n_tt = 1:numel(tt_types)
    if sum(tt_types{n_tt} == '+') || sum(sum(tt_types{n_tt} == ['12']'))==2
        figure;
        sp1 = subplot(2,1,1);
        f_plot_dset_deets(tt_types_data{n_tt}, ops, sp1);
        title(sprintf('data %s, %s',param_name, tt_types{n_tt}), 'interpreter', 'none');
        subplot(2,1,2);
        p_val_mat = f_get_tt_stats(tt_types_data{n_tt});
        imagesc(p_val_mat);
        caxis([0 1]);
    end
end

end