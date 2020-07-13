function f_mpl_spatial_comp(data, ops)
disp('Ensemble analysis...');
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond}; 
    for n_dset = 1:data.(cond_name).num_dsets
        for n_pl = 1:data.(cond_name).num_planes(n_dset)
            A = data.(cond_name).OA_data{n_dset,n_pl}.est.A(:,data.(cond_name).OA_data{n_dset,n_pl}.proc.idx_components);
            A = reshape(A,256,256,[]);
            figure; imagesc(sum(A,3)); title(['Plane ' num2str(n_pl)]);
        end
    end
end
end