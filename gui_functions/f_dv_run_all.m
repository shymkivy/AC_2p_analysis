function f_dv_run_all(app)

num_data = size(app.data,1);
params = f_dv_gather_params(app);

if strcmpi(app.RunallDropDown.Value, 'Stats')
    fprintf('Running all stats, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            params.ddata = ddata;
            if isempty(app.data(n_dset,:).stats{n_pl})
                params.cdata = f_dv_compute_cdata(app, params);
                app.data(n_dset,:).stats{n_pl} = f_dv_compute_stats_core(app, params);
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'dim_red_pca')
    fprintf('Running all dim red pca, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            if isempty(app.data(n_dset,:).data_dim_pca{n_pl})
                params.cdata = f_dv_compute_cdata(app, params);
                app.data(n_dset,:).data_dim_pca{n_pl} = f_dv_estimate_dim_pca_core(params);
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'dim_red_cv')
    fprintf('Running all dim red CV, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            if isempty(app.data(n_dset,:).data_dim_cv{n_pl})
                params.cdata = f_dv_compute_cdata(app, params);
                app.data(n_dset,:).data_dim_cv{n_pl} = f_dv_estimate_dim_cv_core(params);
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'ensemble_extract')
end
disp('Done')
end