% load decoding data from AC_data/decoder folder

trial_data_sort3n = trial_data_sort3 - mean(trial_data_sort3(:));
trial_data_sort3n = trial_data_sort3n/std(trial_data_sort3n(:));

trial_data_sort3nt = trial_data_sort3n';


trial_types2_shuff = trial_types2(randperm(400));



trial_data_sort3nt_past = trial_data_sort3nt(2:end,:);

trial_types2_past = trial_types2(1:end-1);


trial_data_sort3nt_past2 = trial_data_sort3nt(3:end,:);

trial_types2_past2 = trial_types2(1:end-2);





