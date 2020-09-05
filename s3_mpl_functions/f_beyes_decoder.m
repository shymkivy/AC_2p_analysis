function f_beyes_decoder()

tt = [4 5 6 170, 270];

num_cells = size(traces,1);

trials_ind = logical(sum(trial_types==tt,2));

traces2 = traces(:,trials_ind);
trial_types2 = trial_types(trials_ind);




vec_means = zeros(num_cells, numel(tt));
for n_tt = 1:numel(tt)
   tt_ind = trial_types2==tt(n_tt);
   vec_means(:,n_tt) = mean(traces2(:,tt_ind),2);
end

vec_dist = pdist2(vec_means', traces2', 'cosine');

ven_dist_mean = zeros(numel(tt),1);
ven_dist_std = zeros(numel(tt),1);
for n_tt = 1:numel(tt)
    tt_ind = trial_types2==tt(n_tt);
    ven_dist_mean(n_tt) = mean(vec_dist(n_tt,tt_ind));
    ven_dist_std(n_tt) = std(vec_dist(n_tt,tt_ind));
end

num_trials = numel(trial_types2);
vec_decoding = zeros(num_trials,1);
for n_trial = 1:num_trials
    tr_dist = pdist2(vec_means', traces2(:,n_trial)', 'cosine');
    trial_p = exp(-((tr_dist - ven_dist_mean).^2)./(2*ven_dist_std.^2));
    [~, dec_ident] = max(trial_p);
    vec_decoding(n_trial) = tt(dec_ident);
end

sum(trial_types2 == vec_decoding)/num_trials

