function data_out = f_evaluate_ens_result(ens_output, ground_truth)

[gt_mat, gt_list] = if_list_2_mat(ground_truth, []);
[data_ens_mat, ens_list] = if_list_2_mat(ens_output, []);


SI_gt_data = similarity_index(gt_mat, data_ens_mat);

figure;
subplot(2,2,1);
imagesc(SI_gt_data)
ylabel('ground truth ens');
xlabel('data ens');
title('clust pre');

%[~, prem_ind] = max(SI_gt_data);

%SI_gt_data2 = SI_gt_data;
%sizSI = size(SI_gt_data2);
% clust_prem_ind = zeros(numel(gt_list),1);
% for n_tt = 1:(numel(gt_list)-1)
%     [~, temp_ind] = max(SI_gt_data2(:));
%     [gt_ind, dt_ind] = ind2sub(sizSI, temp_ind);
%     clust_prem_ind(gt_ind) = dt_ind;
%     SI_gt_data2(:,dt_ind) = 0;
% end
%clust_prem_ind(clust_prem_ind==0) = find(~sum(ens_list' == clust_prem_ind(clust_prem_ind>0)));

all_perms = perms(1:numel(gt_list));
num_perms = size(all_perms,1);
acc1 = zeros(num_perms,1);
for n_pr = 1:size(all_perms,1)
    temp_SI = similarity_index(gt_mat, data_ens_mat(all_perms(n_pr,:),:));
    acc1(n_pr) = sum(diag(temp_SI));
end
[~, max_perm_ind] = max(acc1);

best_perm = all_perms(max_perm_ind,:);
data_ens_mat_sort = data_ens_mat(best_perm,:);

SI_gt_data3 = similarity_index(gt_mat, data_ens_mat_sort);

subplot(2,2,2);
imagesc(SI_gt_data3);
ylabel('ground truth ens');
xlabel('data ens');
title('clust aligned');

subplot(2,2,3);
bar(diag(SI_gt_data3)); axis tight;
ylim([0 1]);
xlabel('cluster num');
ylabel('overlap accuracy');
title('detection accuracy');

subplot(2,2,4);
bar([sum(gt_mat,2), sum(data_ens_mat_sort,2)]); axis tight;
xlabel('cluster num');
ylabel('num cells');
title('ens size');

data_out.ens_perm_ind = best_perm;
data_out.accuracy = diag(SI_gt_data3);
data_out.ground_truth_cat_order = gt_list;
data_out.data_cat_order = ens_list;

end

function [ens_mat, ens_types] = if_list_2_mat(ens_list, ens_types)

if exist('ens_types', 'var') || isempty(ens_types)
    ens_types = sort(unique(ens_list));
end

ens_mat = zeros(numel(ens_types), numel(ens_list));
for n_cl=1:numel(ens_types)
    ens_mat(n_cl,:) = ens_list == ens_types(n_cl);
end

end