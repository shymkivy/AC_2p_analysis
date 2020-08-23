function f_plot_dset_param_v_ncell(data, param_name, max_cells, ops)

%% collect data
tt_types = unique({data.trial_type_tag})';
tt_types_data = cell(numel(tt_types),1);
for n_tt = 1:numel(tt_types)
    trial_type_tag = tt_types{n_tt};
    tt_cond_data = cell(numel(ops.regions_to_analyze),1);
    dim_est_st2 = data(strcmpi({data.trial_type_tag}, trial_type_tag));
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        dim_data1 = dim_est_st2(strcmpi({dim_est_st2.cond_name}, cond_name));    
        num_dsets = max([dim_data1.n_dset]);
        tt_data = cell(num_dsets,1);
        empty_dset = false(num_dsets,1);
        for n_dset = 1:num_dsets
            dim_data2 = dim_data1([dim_data1.n_dset] == n_dset);
            if numel(dim_data2)
                data1 = [[dim_data2.num_cells_samp]' , [dim_data2.(param_name)]'];
                temp_data = permute(mean(reshape(data1, max([dim_data2.n_rep]),[],2),1), [2 3 1]);
                tt_data{n_dset} = temp_data(temp_data(:,1)<=max_cells,:);
            else
                empty_dset(n_dset) = 1;
            end
        end
        tt_cond_data{n_cond} = tt_data(~empty_dset);
    end
    tt_types_data{n_tt,1} = tt_cond_data;
end

%% combine over dsets
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

%% plot stuff

for n_tt = 1:numel(tt_types)
    if sum(tt_types{n_tt} == '+') || sum(sum(tt_types{n_tt} == ['12']'))==2
        f1 = figure; hold on;
        pl1 = cell(numel(ops.regions_to_analyze),1);
        dset_data_mean = cell(numel(ops.regions_to_analyze),1);
        for n_cond = 1:numel(ops.regions_to_analyze)
            num_dsets = numel(tt_types_data{n_tt}{n_cond});
            max_size = 0;
            for n_pl = 1:num_dsets
                max_size = max(max_size, size(tt_types_data{n_tt}{n_cond}{n_pl},1));
            end
            dset_data = nan(max_size,2,num_dsets);
            for n_dset = 1:num_dsets
                temp_data = tt_types_data{n_tt}{n_cond}{n_dset};
                dset_data(1:size(temp_data,1),:,n_dset) = temp_data;
                figure(f1);
                p1 = plot(temp_data(:,1),temp_data(:,2), 'Color', ops.cond_colors{n_cond}, 'LineWidth', 0.1);
                p1.Color(4) = 0.5;
            end
            dset_data_mean{n_cond} = nanmean(dset_data,3);
        end
        for n_cond = 1:numel(ops.regions_to_analyze)
            figure(f1);
            pl1{n_cond} = plot(dset_data_mean{n_cond}(:,1),dset_data_mean{n_cond}(:,2), 'Color', ops.cond_colors{n_cond}, 'LineWidth', 3);
        end
        figure(f1);
        xlabel('Number of cells');
        ylabel('Estimated rank');
        %xlim([0 300]);
        %ylim([0 25]);
        legend([pl1{1} pl1{2} pl1{3} pl1{4}], ops.regions_to_analyze,'Location','southeast');
        title(sprintf('%s, %s', tt_types{n_tt}, param_name), 'Interpreter', 'none');
    end
end




end