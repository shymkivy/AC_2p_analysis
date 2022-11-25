function [trial_data_sort_wctx, trial_types_wctx, trial_types_idx_wctx] = f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, ops)

trial_types_idx_wctx = (1:numel(trial_types))';
trial_types_wctx = trial_types;
if ~isempty(trial_data_sort)
    trial_data_sort_wctx = trial_data_sort;
else
    trial_data_sort_wctx = [];
end

if ~isempty(mmn_freq)
    % add cont 18, 150 (mmn_freq(2))
    cont_ind = logical(trial_types == mmn_freq(2));
    cont_ind_list = find(cont_ind);
    % add redf 19, 260 (mmn_freq(2))
    redf_ind = logical(sum(trial_types == (200 + ops.redundent_pool_trials),2));
    redf_ind_list = find(redf_ind);
    % add contf 28 (mmn_freq(1))
    cont2_ind = logical(trial_types == mmn_freq(1));
    cont2_ind_list = find(cont2_ind);
    % add red 29   (mmn_freq(1))
    red_ind = logical(sum(trial_types == (100 + ops.redundent_pool_trials),2));
    red_ind_list = find(red_ind);

    trial_types_wctx = cat(1, trial_types_wctx, 150*ones(sum(cont_ind),1));
    trial_types_wctx = cat(1, trial_types_wctx, 260*ones(sum(redf_ind),1));
    trial_types_wctx = cat(1, trial_types_wctx, 250*ones(sum(cont2_ind),1));
    trial_types_wctx = cat(1, trial_types_wctx, 160*ones(sum(red_ind),1));

    trial_types_idx_wctx = cat(1, trial_types_idx_wctx, cont_ind_list);
    trial_types_idx_wctx = cat(1, trial_types_idx_wctx, redf_ind_list);
    trial_types_idx_wctx = cat(1, trial_types_idx_wctx, cont2_ind_list);
    trial_types_idx_wctx = cat(1, trial_types_idx_wctx, red_ind_list); 
    
    if ~isempty(trial_data_sort)
        trial_data_sort_wctx = cat(3, trial_data_sort_wctx, trial_data_sort(:,:,cont_ind));
        trial_data_sort_wctx = cat(3, trial_data_sort_wctx, trial_data_sort(:,:,redf_ind));
        trial_data_sort_wctx = cat(3, trial_data_sort_wctx, trial_data_sort(:,:,cont2_ind));
        trial_data_sort_wctx = cat(3, trial_data_sort_wctx, trial_data_sort(:,:,red_ind));
    end
end
end