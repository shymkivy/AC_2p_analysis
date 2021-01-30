function f_plot_cond_reliability(data, ops)

%e_colors = {'b', 'r', 'g', 'k'};

for n_tt = 1:numel(ops.dred_params.trial_types_to_dred)
    figure; hold on;
    full_list = cell(numel(ops.regions_to_analyze),1);
    all_trials = cell(numel(ops.regions_to_analyze),1);
    num_cell_cond = zeros(numel(ops.regions_to_analyze),1);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        cdata = data(strcmpi(data.area, cond_name),:);

        dr_params.cond_name = cond_name;
        dr_params.colors_clust = cat(2,ops.colors_list,ops.colors_list,ops.colors_list);
        dr_params.tt_to_dred_input = ops.dred_params.trial_types_to_dred{n_tt};
        dr_params.num_clust = ops.dred_params.hclust.num_clust{n_tt};
        
        reliab_list = cell(numel(cdata.area),1);
        for n_dset = 1:numel(cdata.area)
            [tn_to_dred, trial_type_tag] = f_select_trial_type(dr_params.tt_to_dred_input, cdata, n_dset, ops);
            
            reliab_list_temp = cell(numel(tn_to_dred),1);
            for n_tt2 = 1:numel(tn_to_dred)
                resp_cells = logical(cdata.peak_tuned_trials_full{n_dset}(:,tn_to_dred(n_tt2)));
                cell_reliab = cdata.peak_tuned_trials_full_reliab{n_dset}(:,tn_to_dred(n_tt2));
                reliab_list_temp{n_tt2} = cell_reliab(resp_cells);
            end
            
            reliab_list{n_dset} = cat(1,reliab_list_temp{:});
        end
        reliab_list = cat(1,reliab_list{:});
        
        num_cell_cond(n_cond) = numel(reliab_list);
         
        full_list{n_cond} = reliab_list;
        %[f, x] = ksdensity(dd_resp_reliab, 'Bandwidth',.07);

        %histogram(dd_resp_reliab, 10, 'BinLimits',[0,1], 'FaceColor','none','EdgeColor',e_colors{n_cond},'Normalization','probability', 'LineWidth', 2,'DisplayStyle','stairs')
        
        x = sort(reliab_list);
        f = (1:numel(reliab_list))/numel(reliab_list);

        
        %[f, x] = ecdf(reliab_list);
        plot(x,f, 'color', ops.cond_colors{n_cond}, 'LineWidth', 2);

        all_trials{n_cond} = cat(1,cdata.peak_tuned_trials_full_reliab{:});
    end
    %% make shuff list
    
    full_list = cat(1,full_list{:});
    
    num_repeats = 1000;
    num_samp = round(mean(num_cell_cond));
    
    samp_list = if_sample_trials(full_list, num_repeats, num_samp);
    
    shadedErrorBar(mean(samp_list,1),(1:num_samp)/num_samp, std(samp_list, [], 1))
    
    all_trials = cat(1,all_trials{:});
    samp_list2 = if_sample_trials(all_trials(:), num_repeats, num_samp);
%   shadedErrorBar(mean(samp_list2,1),(1:num_samp)/num_samp, std(samp_list, [], 1))
    
    %plot(x,f, 'LineStyle', '--', 'LineWidth', 2, 'Color', [.8 .8 .8]);
    
    axis tight;
    legend([ops.regions_to_analyze {'Shuff'}], 'Location', 'southeast');
    xlabel('response reliability');
    ylabel('fraction');
    title(sprintf('Cell reliability ecdf, trials %s', trial_type_tag));
end


end

function samp_list = if_sample_trials(full_list, num_repeats, num_samp)

samp_list = zeros(num_repeats, num_samp);
for n_rand = 1:num_repeats
    samp_list(n_rand,:) = sort(randsample(full_list, num_samp));
end

end