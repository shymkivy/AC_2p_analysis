function f_plot_stim_vec_dist_v_ncell(trial_ave_vec, tt_ind, stim_ind, ops, max_cells)
start_range = 1;
cell_int = 2;
num_samp = 10;
max_cells = 100;

max_len = 0;
corr_list = cell(numel(ops.regions_to_analyze),1);
for n_cond = 1:numel(ops.regions_to_analyze)
    num_dsets = size(trial_ave_vec{n_cond},2);
    corr_list{n_cond} = zeros(num_dsets,1);
    corr_list2 = cell(size(stim_ind,1), num_dsets);
    for n_stim = 1:size(stim_ind,1)
        
        for n_dset = 1:num_dsets
            temp_data = trial_ave_vec{n_cond}{tt_ind,n_dset};
            num_cells = min(size(temp_data,1),max_cells);
            cell_range = start_range:cell_int:num_cells;
            num_range = numel(cell_range);
            max_len = max(max_len, num_range);
            corr_list3 = zeros(num_range, num_samp);
            for n_range = 1:num_range
                for n_rep = 1:num_samp
                    ls1 = randsample(num_cells,cell_range(n_range));
                    
                    corr_list3(n_range,n_rep) = 1-pdist2(temp_data(ls1, stim_ind(n_stim,1))',temp_data(ls1, stim_ind(n_stim,2))','cosine');
                end
            end
            corr_list2{n_stim, n_dset} = nanmean(corr_list3,2);
        end
        
    end
    corr_list{n_cond} = corr_list2(:);
end


%% plot stuff


f1 = figure; hold on;
pl1 = cell(numel(ops.regions_to_analyze),1);
dset_data_mean = cell(numel(ops.regions_to_analyze),1);
for n_cond = 1:numel(ops.regions_to_analyze)
    num_dsets = numel(corr_list{n_cond});
    dset_data = nan(max_len,2,num_dsets);
    for n_dset = 1:num_dsets
        temp_data = corr_list{n_cond}{n_dset};
        dset_data(1:numel(temp_data),1,n_dset) = 1:cell_int:((numel(temp_data))*cell_int);
        dset_data(1:numel(temp_data),2,n_dset) = temp_data;
        figure(f1);
        p1 = plot(dset_data(1:numel(temp_data),1,n_dset),dset_data(1:numel(temp_data),2,n_dset), 'Color', ops.cond_colors{n_cond}, 'LineWidth', 0.1);
        p1.Color(4) = 0.3;
    end
    dset_data_mean{n_cond} = nanmean(dset_data,3);
end
for n_cond = 1:numel(ops.regions_to_analyze)
    figure(f1);
    
    pl1{n_cond} = plot(dset_data_mean{n_cond}(:,1),dset_data_mean{n_cond}(:,2), 'Color', ops.cond_colors{n_cond}, 'LineWidth', 3);
end
figure(f1); axis tight;
xlabel('Number of cells');
ylabel('Estimated similarity');
%xlim([0 300]);
%ylim([0 25]);
legend([pl1{1} pl1{2} pl1{3} pl1{4}], ops.regions_to_analyze,'Location','southeast');



end