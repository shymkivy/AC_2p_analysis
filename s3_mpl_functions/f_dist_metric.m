function f_dist_metric(trial_peaks, trial_types, dr_params)

tt_to_dred = dr_params.tt_to_dred;

trial_peaks2 = trial_peaks(dr_params.hclust_out_cell.dend_order,:);
num_cells = size(trial_peaks2,1);

vec1 = trial_peaks2(:,trial_types == tt_to_dred(1));
vec2 = trial_peaks2(:,trial_types == tt_to_dred(2));

figure;
subplot(1,4,1);
imagesc(vec1);
subplot(1,4,2)
plot(mean(vec1,2),1:num_cells)
subplot(1,4,3)
imagesc(vec2)
subplot(1,4,4)
plot(mean(vec2,2),1:num_cells)

dist1 = 1 - pdist2(vec1', vec2', 'cosine');

f_ensemble_analysis_peaks(trial_peaks_dred_sort(:,tt_index1), trial_types_dred_sort(tt_index1), dr_params, ops);

end