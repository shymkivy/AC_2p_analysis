function align_mat_out = f_align_vec(vec_templ, vec_al, metric)
% input vectors (num_cluster x feature(axes))

%vec_al_norm = vec_al./vecnorm(vec_al')';

SI = 1 - pdist2(vec_templ, vec_al, metric);

% figure; imagesc(SI)
% ylabel('Vec 1');
% xlabel('Vec 2');
% title('pre')

v1_d = size(vec_templ,1);
v2_d = size(vec_al,1);

[max_size, max_ind] = max([v1_d, v2_d]);
min_size = min([v1_d, v2_d]);
all_perms = perms(1:max_size);
num_perms = size(all_perms,1);
acc1 = zeros(num_perms,1);
for n_pr = 1:size(all_perms,1)
    if max_ind == 2
        temp_SI = SI(:,all_perms(n_pr,:));
    else
        temp_SI = SI(all_perms(n_pr,:),:);
    end
    if min_size == 1
        acc1(n_pr) = temp_SI(1,1);
    else
        acc1(n_pr) = sum(diag(temp_SI));
    end
    
end
[~, max_perm_ind] = max(acc1);

best_perm = all_perms(max_perm_ind,:);

if max_ind == 2
    best_SI = SI(:,best_perm);
else
    best_SI = SI(best_perm,:);
end

% figure; imagesc(best_SI)
% ylabel('Vec 1');
% xlabel('Vec 2');
% title('post');

align_mat1 = (1:v1_d)';
align_mat2 = (1:v2_d)';

if max_ind == 2
    align_mat2 = align_mat2(best_perm);
else
    align_mat1 = align_mat1(best_perm);
end

align_mat_out = [align_mat1(1:min_size), align_mat2(1:min_size)];

[~, B] = sort(align_mat_out(:,2));
align_mat_out = align_mat_out(B,:);

end