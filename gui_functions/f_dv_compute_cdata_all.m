function f_dv_compute_cdata_all(app)

num_dsets = size(app.data,1);
if ~sum(strcmpi(app.data.Properties.VariableNames, 'cdata'))
    max_planes = max(app.data.num_planes);
    app.data.cdata = cell(num_dsets,max_planes);
end

wb = f_waitbar_initialize(app, 'Computing cdata...');
params = f_dv_gather_params(app);
for n_dset = 1:num_dsets
    f_waitbar_update(wb, n_dset/num_dsets, sprintf('Computing cdata %d/%d...', n_dset, num_dsets));
    params.n_dset = n_dset;
    cdata = cell(1, app.data(n_dset,:).num_planes);
    for n_pl = 1:app.data(n_dset,:).num_planes
        params.n_pl = n_pl;
        cdata{n_pl} = f_dv_compute_cdata(app, params);
    end
    app.data(n_dset,:).cdata = cdata;  
end
f_waitbar_close(wb);

end