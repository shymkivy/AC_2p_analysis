function f_dv_run_all_stats(app)

num_data = size(app.data,1);

params = f_dv_gather_params(app);

fprintf('Dset_/%d: ', num_data)
for n_dset = 1:num_data
    fprintf('%d..', n_dset);
    params.n_dset = n_dset;
    ddata = app.data(n_dset,:);
    if isempty(app.data(n_dset,:).stats{1})
        app.data(n_dset,:).stats{1} = cell(ddata.num_planes,1);
    end 
    for n_pl = 1:ddata.num_planes
        params.n_pl = n_pl;
        params.cdata = f_dv_compute_cdata(app, params);
        if isempty(app.data(n_dset,:).stats{1}{n_pl})
            app.data(n_dset,:).stats{1}{n_pl} = f_dv_compute_stats(app,params);
        end
    end
end
fprintf('\n');

end