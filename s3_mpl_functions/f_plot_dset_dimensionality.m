function f_plot_dset_dimensionality(data, ops)

%% total
tt_types = fields(data);
tt_types_data = cell(numel(tt_types),2);
for n_tt = 1:numel(tt_types)
    tt_dim_cond = cell(numel(ops.regions_to_analyze),1);
    tt_dim_cond_shuff = cell(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        chdata = data.(tt_types{n_tt}).(cond_name);
        
        num_dsets = numel(chdata.hclust_out_tr);
        tt_dim = zeros(num_dsets,1);
        tt_dim_shuff = zeros(num_dsets,1);
        for n_dset = 1:num_dsets
            tt_dim(n_dset) = chdata.data_dim_est_full{n_dset}.dimensionality_total;
            tt_dim_shuff(n_dset) = chdata.data_dim_est_full{n_dset}.dimensionality_total_shuff;
        end
        tt_dim_cond{n_cond} = tt_dim;
        tt_dim_cond_shuff{n_cond} = tt_dim_shuff;
    end
    tt_types_data{n_tt,1} = tt_dim_cond;
    tt_types_data{n_tt,2} = tt_dim_cond_shuff;
end
if sum(strcmpi(tt_types, 'dd1')) && sum(strcmpi(tt_types, 'dd2'))
    size1 = numel(tt_types)+1;
    for n_cond = 1:numel(ops.regions_to_analyze)
        temp1 = [tt_types_data{strcmpi(tt_types, 'dd1'),1}{n_cond}; tt_types_data{strcmpi(tt_types, 'dd2'),1}{n_cond}];
        %temp2 = [tt_types_data{strcmpi(tt_types, 'dd1'),2}{n_cond}; tt_types_data{strcmpi(tt_types, 'dd2'),2}{n_cond}];
        tt_types_data{size1,1}{n_cond,1} = temp1;
        %tt_types_data{size1,2}{n_cond,1} = temp2;
    end 
    tt_types = cat(1, tt_types, 'dd1+2');
end
for n_tt = 1:numel(tt_types)
    f_plot_dset_deets(tt_types_data{n_tt}, ops);
    %f_plot_dset_deets(cat(2,tt_types_data{n_tt,:}), ops);
    title(sprintf('data total dimensionality, %s', tt_types{n_tt}));
end

%% correlations
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

%%
tt_types = fields(data);
tt_types_data = cell(numel(tt_types),2);
for n_tt = 1:numel(tt_types)
    tt_dim_cond_shuff = cell(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        chdata = data.(tt_types{n_tt}).(cond_name);
        
        num_dsets = numel(chdata.hclust_out_tr);
        tt_dim = zeros(num_dsets,1);
        tt_dim_shuff = zeros(num_dsets,1);
        for n_dset = 1:num_dsets
            tt_dim_shuff(n_dset) = chdata.data_dim_est_full{n_dset}.dimensionality_total_shuff;
        end
        tt_dim_cond{n_cond} = tt_dim;
        tt_dim_cond_shuff{n_cond} = tt_dim_shuff;
    end
    tt_types_data{n_tt,1} = tt_dim_cond_shuff;
end
if sum(strcmpi(tt_types, 'dd1')) && sum(strcmpi(tt_types, 'dd2'))
    size1 = numel(tt_types)+1;
    for n_cond = 1:numel(ops.regions_to_analyze)
        temp1 = [tt_types_data{strcmpi(tt_types, 'dd1'),1}{n_cond}; tt_types_data{strcmpi(tt_types, 'dd2'),1}{n_cond}];
        %temp2 = [tt_types_data{strcmpi(tt_types, 'dd1'),2}{n_cond}; tt_types_data{strcmpi(tt_types, 'dd2'),2}{n_cond}];
        tt_types_data{size1,1}{n_cond,1} = temp1;
        %tt_types_data{size1,2}{n_cond,1} = temp2;
    end 
    tt_types = cat(1, tt_types, 'dd1+2');
end
for n_tt = 1:numel(tt_types)
    f_plot_dset_deets(tt_types_data{n_tt}, ops);
    title(sprintf('data total dimensionality shuff, %s', tt_types{n_tt}));
end

end