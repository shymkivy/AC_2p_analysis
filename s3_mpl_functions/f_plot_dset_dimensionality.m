function f_plot_dset_dimensionality(data, ops)

tt_types = fields(data);
tt_types_data = cell(numel(tt_types),1);
for n_tt = 1:numel(tt_types)
    tt_dim_cond = cell(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        chdata = data.(tt_types{n_tt}).(cond_name);
        
        num_dsets = numel(chdata.hclust_out_tr);
        tt_dim = zeros(num_dsets,1);
        for n_dset = 1:num_dsets
            tt_dim(n_dset) = chdata.data_dim_est_full{n_dset}.dimensionality_total;
        end
        tt_dim_cond{n_cond} = tt_dim;
    end
    tt_types_data{n_tt} = tt_dim_cond;
end
if sum(strcmpi(tt_types, 'dd1')) && sum(strcmpi(tt_types, 'dd2'))
    size1 = numel(tt_types_data)+1;
    for n_cond = 1:numel(ops.regions_to_analyze)
        tt_types_data{size1,1}{n_cond,1} = [tt_types_data{strcmpi(tt_types, 'dd1')}{n_cond}; tt_types_data{strcmpi(tt_types, 'dd2')}{n_cond}];
    end 
    tt_types = cat(1, tt_types, 'dd1+2');
end

for n_tt = 1:numel(tt_types)
    f_plot_dset_deets(tt_types_data{n_tt}, ops);
    title(sprintf('data total dimensionality, %s', tt_types{n_tt}));
end

tt_types = fields(data);
tt_types_data = cell(numel(tt_types),1);
for n_tt = 1:numel(tt_types)
    tt_dim_cond = cell(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        chdata = data.(tt_types{n_tt}).(cond_name);
        
        num_dsets = numel(chdata.hclust_out_tr);
        tt_dim = zeros(num_dsets,1);
        for n_dset = 1:num_dsets
            tt_dim(n_dset) = chdata.data_dim_est_full{n_dset}.dimensionality_corr;
        end
        tt_dim_cond{n_cond} = tt_dim;
    end
    tt_types_data{n_tt} = tt_dim_cond;
end
if sum(strcmpi(tt_types, 'dd1')) && sum(strcmpi(tt_types, 'dd2'))
    size1 = numel(tt_types_data)+1;
    for n_cond = 1:numel(ops.regions_to_analyze)
        tt_types_data{size1,1}{n_cond,1} = [tt_types_data{strcmpi(tt_types, 'dd1')}{n_cond}; tt_types_data{strcmpi(tt_types, 'dd2')}{n_cond}];
    end 
    tt_types = cat(1, tt_types, 'dd1+2');
end

for n_tt = 1:numel(tt_types)
    f_plot_dset_deets(tt_types_data{n_tt}, ops);
    title(sprintf('data correlations dimensionality, %s', tt_types{n_tt}));
end


end