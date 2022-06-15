function f_dv_comp_tuning_similarity(app)

ddata = app.ddata;
num_planes = ddata.num_planes;

idx1 = logical(strcmpi(ddata.mouse_id, app.data.mouse_id).*(ddata.FOV_num == app.data.FOV_num));
data2 = app.data(idx1,:);

num_dsets = sum(idx1);

match_stats_id = data2.registration{1}.match_stats_id;

trials1 = 1:10;

vec_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    vec_all2 = cell(num_planes,1);
    for n_pl = 1:num_planes
        stats1 = data2.stats{n_dset,n_pl};
        match_stats_id2 = match_stats_id{n_pl};
        
        peak_vals = stats1.peak_val_all(match_stats_id2(:,n_dset), trials1);
        
        vec_all2{n_pl} = peak_vals;
    end
    vec_all{n_dset} = cat(1, vec_all2{:});
end

vec_all3 = reshape(cat(3, vec_all{:}), [], num_dsets);

si_cos = 1 - pdist2(vec_all3', vec_all3', 'cosine');

si_corr = 1 - pdist2(vec_all3', vec_all3', 'correlation');

figure; imagesc(si_cos); axis equal tight;
title(sprintf('%s; population trial ave cosine similarity', ddata.mouse_id{1}), 'interpreter', 'none');

figure; imagesc(si_corr); axis equal tight;
title(sprintf('%s; population trial ave correlation', ddata.mouse_id{1}), 'interpreter', 'none');

end