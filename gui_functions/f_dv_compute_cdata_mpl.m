function cdata = f_dv_compute_cdata_mpl(ddata, params)

cdata = cell(1, ddata.num_planes);
for n_pl = 1:ddata.num_planes
    fprintf('%s, pl%d\n', ddata.dset_name_full{1}, n_pl);
    params.n_pl = n_pl;
    cdata{n_pl} = f_dv_compute_cdata(ddata, params);
end

end