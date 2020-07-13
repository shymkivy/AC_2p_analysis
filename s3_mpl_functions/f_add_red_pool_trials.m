function [trial_data_sort_poolred, trial_types_poolred] = f_add_red_pool_trials(trial_data_sort, trial_types, ops)

redf_ind = logical(sum(trial_types == (200 + ops.redundent_pool_trials),2));
trial_data_sort_poolred = cat(3, trial_data_sort, trial_data_sort(:,:,redf_ind));
trial_types_poolred = cat(1, trial_types, 260*ones(sum(redf_ind),1));

red_ind = logical(sum(trial_types == (100 + ops.redundent_pool_trials),2));
trial_data_sort_poolred = cat(3, trial_data_sort_poolred, trial_data_sort(:,:,red_ind));
trial_types_poolred = cat(1, trial_types_poolred, 160*ones(sum(red_ind),1));

end