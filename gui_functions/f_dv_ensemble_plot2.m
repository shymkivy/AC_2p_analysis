function f_dv_ensemble_plot2(app)

1

ddata = app.ddata;
cdata = f_dv_get_cdata(app);

ens = ddata.ensembles{1};
ens_stats = ddata.ensemble_stats{1};


A_all = cell(ddata.num_planes, 1);
cell_plane_idx = cell(ddata.num_planes, 1);
for n_pl = 1:ddata.num_planes
    est = ddata.OA_data{n_pl}.est;
    proc = ddata.OA_data{n_pl}.proc;
    cell_plane_idx{n_pl} = ones(sum(proc.comp_accepted),1)*n_pl;
    A_all{n_pl} = est.A(:,proc.comp_accepted);
end

cell_plane_idx = cat(1,cell_plane_idx{:});


planes_A = cell(ddata.num_planes,1);

    
    
for n_pl = 1:ddata.num_planes
    for n_ens = 1:ens.num_ens_comps
        
        coeffs1 = ens.coeffs(:,n_ens);
        thresh = 2*std(coeffs1);
        
        x1 = coeffs1(cell_plane_idx==n_pl);
        idx1 = coeffs1(cell_plane_idx==n_pl)>thresh;
        im1 = repmat(full(A_all{n_pl}),1,1,3);
        im1(:,idx1,:) = im1(:,idx1,:).*repmat(x1(idx1)',1,1,3).*reshape([1 0 0],1,1,3)/10;
        planes_A{n_pl} = squeeze(sum(reshape(im1, est.dims(1), est.dims(2), [], 3),3));
    end

end

planes_A_all = cat(2, planes_A{:});


figure; imagesc(planes_A_all*3)

reshape(sum(A_all{1},3), 256, 256)


sum(coeffs1>thresh)

sort(coeffs1)


end