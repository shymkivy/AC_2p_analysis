function f_plot_dim_est_data(dim_est_st, ops)

for n_tt = 1:numel(ops.dred_params.trial_types_to_dred)
    if isnumeric(ops.dred_params.trial_types_to_dred{n_tt})
        trial_type_tag = ['trials_' trial_types_to_dred];
    else
        trial_type_tag = ops.dred_params.trial_types_to_dred{n_tt};
    end
    
    dim_est_st2 = dim_est_st(strcmpi({dim_est_st.trial_type_tag}, trial_type_tag));
    colors1 = ops.cond_colors;%{[.5 .5 1], [1 .5 .5], [.5 1 .5]};
    colors2 = ops.cond_colors;%{[0 0 1], [1 0 0], [0 1 0], [1 1 ]};
    figure; 
    s1 = subplot(1,2,1); hold on;
    s2 = subplot(1,2,2); hold on;
    pl1 = cell(numel(ops.regions_to_analyze),1);
    pl2 = cell(numel(ops.regions_to_analyze),1);
    dset_data_mean = cell(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        dim_data1 = dim_est_st2(strcmpi({dim_est_st2.cond_name}, cond_name));
        for ii = 1:numel(dim_data1)
            dim_data1(ii).d_explained = sum(dim_data1(ii).d_explained);
        end
        dset_data_size = 0;
        for n_dset = 1:max([dim_data1.n_dset])
            dim_data2 = dim_data1([dim_data1.n_dset] == n_dset);
            data1 = [[dim_data2.num_cells_samp]' , [dim_data2.dimensionality_est]', [dim_data2.d_explained]'];
            data2 = permute(mean(reshape(data1, max([dim_data2.n_rep]),[],3),1), [2 3 1]);
            dset_data_size = max([dset_data_size, size(data2,1)]);
        end
        dset_data = nan(dset_data_size-1,2,max([dim_data1.n_dset]));
        dset_data_lamb = nan(dset_data_size-1,2,max([dim_data1.n_dset]));
        for n_dset = 1:max([dim_data1.n_dset])
            dim_data2 = dim_data1([dim_data1.n_dset] == n_dset);
            data1 = [[dim_data2.num_cells_samp]' , [dim_data2.dimensionality_est]', [dim_data2.d_explained]'];
            data2 = permute(mean(reshape(data1, max([dim_data2.n_rep]),[],3),1), [2 3 1]);
            subplot(s1)
            %plot(data2(:,1),data2(:,2), 'Color', colors1{n_cond}, 'LineWidth', 0.1);
            subplot(s2)
            %plot(data2(:,2),data2(:,3), 'Color', colors1{n_cond}, 'LineWidth', 0.1);
            dset_data(1:(size(data2,1)-1),:,n_dset) = data2(1:(size(data2,1)-1),1:2);
            dset_data_lamb(1:(size(data2,1)-1),:,n_dset) = data2(1:(size(data2,1)-1),[2, 3]);
        end
        dset_data_mean{n_cond} = nanmean(dset_data,3);
        dset_data_lamb_mean{n_cond} = nanmean(dset_data_lamb,3);
    end
    for n_cond = 1:numel(ops.regions_to_analyze)
        subplot(s1)
        pl1{n_cond} = plot(dset_data_mean{n_cond}(:,1),dset_data_mean{n_cond}(:,2), 'Color', colors2{n_cond}, 'LineWidth', 3);
        subplot(s2)
        pl2{n_cond} = plot(dset_data_lamb_mean{n_cond}(:,1),dset_data_lamb_mean{n_cond}(:,2), 'Color', colors2{n_cond}, 'LineWidth', 3);
    end
    subplot(s1)
    xlabel('Number of cells');
    ylabel('Estimated rank');
    %xlim([0 300]);
    %ylim([0 25]);
    legend([pl1{1} pl1{2} pl1{3} pl1{4}], ops.regions_to_analyze,'Location','southeast');
    title('Dimensionality of correlations', 'Interpreter', 'none');
    subplot(s2)
    xlabel('Data Rank');
    ylabel('Variance explained (%)');
    legend([pl2{1} pl2{2} pl2{3} pl1{4}], ops.regions_to_analyze,'Location','southeast');
    title('Variance in low rank data');
    suptitle(['trials analyzed ' trial_type_tag])
end

end