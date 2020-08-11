function f_measure_temp_similarity(scores, trial_types, tt_to_dred)

num_rep = 10000;
figure; hold on;
for n_tt = 1:numel(tt_to_dred)
    tt_list = find(trial_types == tt_to_dred(n_tt));
    dist1 = zeros((numel(tt_list)-1),1);
    for n_rep = 1:(numel(tt_list)-1)
        dist1(n_rep) = (scores(:,tt_list(n_rep))/norm(scores(:,tt_list(n_rep))))'*(scores(:,tt_list(n_rep+1))/norm(scores(:,tt_list(n_rep+1))));
    end
    [f, x] = ecdf(dist1);
    plot(x, f, 'color', ops.colors_list{n_tt}, 'lineWidth', 2);
end

for n_tt = 1:numel(tt_to_dred)
    dist1 = zeros(num_rep,3);
    tt_list = find(trial_types == tt_to_dred(n_tt));
    for n_rep = 1:num_rep
        samp_tr = randsample(tt_list,2);
        dist1(n_rep,n_tt) = (scores(:,samp_tr(1))/norm(scores(:,samp_tr(1))))'*(scores(:,samp_tr(2))/norm(scores(:,samp_tr(2))));
    end
    [f, x] = ecdf(dist1(:,n_tt));
    plot(x, f, 'color', ops.colors_list{n_tt}, 'LineStyle', '--', 'lineWidth', 2);
end
title(sprintf('%s, dset%d, Onset population trajectory structure ECDF',params.cond_name, params.n_dset));
legend('cont', 'red', 'dev', 'rand', 'Location', 'northwest');
xlabel('cosine similarity betweeen trials');
ylabel('Fraction');


end