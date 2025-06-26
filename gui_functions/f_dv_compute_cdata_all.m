function data = f_dv_compute_cdata_all(data, ops)
params = ops.params;

num_dsets = size(data,1);
if ~sum(strcmpi(data.Properties.VariableNames, 'cdata'))
    max_planes = max(data.num_planes);
    data.cdata = cell(num_dsets,max_planes);
end

if ops.waitbar
    if isfield(ops, 'app')
        app = ops.app;
    else
        app = [];
    end
    wb = f_waitbar_initialize(app, 'Computing cdata...');
end
%params = f_dv_gather_params(app);
for n_dset = 1:num_dsets
    if ops.waitbar
        f_waitbar_update(wb, n_dset/num_dsets, sprintf('Computing cdata %d/%d...', n_dset, num_dsets));
    end
    ddata = data(n_dset,:);
    params.n_dset = n_dset;
    data(n_dset,:).cdata = f_dv_compute_cdata_mpl(ddata, params);
end
if ops.waitbar
    f_waitbar_close(wb);
end

end