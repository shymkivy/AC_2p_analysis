function [trial_data_sort_wctx, trial_types_wctx] = f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, ops)

% add cont 18, 150
cont_ind = logical(trial_types == mmn_freq(2));
% add redf 19, 260
redf_ind = logical(sum(trial_types == (200 + ops.redundent_pool_trials),2));
% add contf 28
cont2_ind = logical(trial_types == mmn_freq(1));
% add red 29
red_ind = logical(sum(trial_types == (100 + ops.redundent_pool_trials),2));


trial_types_wctx = cat(1, trial_types, 150*ones(sum(cont_ind),1));
trial_types_wctx = cat(1, trial_types_wctx, 260*ones(sum(redf_ind),1));
trial_types_wctx = cat(1, trial_types_wctx, 250*ones(sum(cont2_ind),1));
trial_types_wctx = cat(1, trial_types_wctx, 160*ones(sum(red_ind),1));

if ~isempty(trial_data_sort)
    trial_data_sort_wctx = cat(3, trial_data_sort, trial_data_sort(:,:,cont_ind));
    trial_data_sort_wctx = cat(3, trial_data_sort_wctx, trial_data_sort(:,:,redf_ind));
    trial_data_sort_wctx = cat(3, trial_data_sort_wctx, trial_data_sort(:,:,cont2_ind));
    trial_data_sort_wctx = cat(3, trial_data_sort_wctx, trial_data_sort(:,:,red_ind));
else
    trial_data_sort_wctx = [];
end

end