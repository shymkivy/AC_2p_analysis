function f_dv_comp_tuning_similarity(app)

ddata = app.ddata;
num_planes = ddata.num_planes;
params = f_dv_gather_params(app);

idx1 = logical(strcmpi(ddata.mouse_id, app.data.mouse_id).*(ddata.FOV_num == app.data.FOV_num));
data2 = app.data(idx1,:);

num_dsets = sum(idx1);

trials1 = 1:10;

[sort_names, sort_idx] = sort(data2.dset_name);

% find index of common cells in all dsets and also any resp cells (or)
pl_cell_idx = cell(num_planes,1);
resp_cells_idx = cell(num_planes,1);
for n_pl = 1:num_planes
    for n_dset = 1:num_dsets
        if n_dset == 1
            num_cells = numel(data2.register_roi{n_dset,n_pl}.fov_cell_idx);
            pl_cell_idx{n_pl} = true(num_cells,1);
            resp_cells_idx{n_pl} = false(num_cells,1);
        end
        reg_cell_idx = ~isnan(data2.register_roi{n_dset,n_pl}.reg_cell_idx);
        pl_cell_idx{n_pl} = and(pl_cell_idx{n_pl}, reg_cell_idx);
        
        stats1 = data2.stats{n_dset,n_pl};
        [resp_cells, ~, ~] = f_dv_get_resp_vals_cells(stats1, trials1, params);
        resp_cells2 = logical(sum(resp_cells,2));
        
        resp_cells_idx{n_pl}(reg_cell_idx) = or(resp_cells_idx{n_pl}(reg_cell_idx), resp_cells2);
    end
end

% gather data
vec_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    vec_all2 = cell(num_planes,1);
    
    for n_pl = 1:num_planes
        stats1 = data2.stats{n_dset,n_pl};
    
        [~, ~, resp_vals] = f_dv_get_resp_vals_cells(stats1, trials1, params);
        
        resp_idx = and(resp_cells_idx{n_pl}, pl_cell_idx{n_pl});
        resp_cells3 = data2.register_roi{n_dset,n_pl}.reg_cell_idx(resp_idx);
        
        resp_vals2 = resp_vals(resp_cells3,:);
        
        vec_all2{n_pl} = resp_vals2;
    end
    vec_all{n_dset} = cat(1, vec_all2{:});
end

vec_all3 = reshape(cat(3, vec_all{:}), [], num_dsets);

vec_all4 = vec_all3(:,sort_idx);

figure; plot(vec_all4)
legend(sort_names)

si_cos = 1 - pdist2(vec_all4', vec_all4', 'cosine');

si_corr = 1 - pdist2(vec_all4', vec_all4', 'correlation');

figure; imagesc(si_cos); axis equal tight;
title(sprintf('%s; population trial ave cosine similarity', ddata.mouse_id{1}), 'interpreter', 'none');

figure; imagesc(si_corr); axis equal tight;
title(sprintf('%s; population trial ave correlation', ddata.mouse_id{1}), 'interpreter', 'none');

end