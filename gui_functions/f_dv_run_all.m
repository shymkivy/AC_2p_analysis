function f_dv_run_all(app)

num_data = size(app.data,1);
params = f_dv_gather_params(app);

% reset in case
%max_planes = max(app.data.num_planes);
%app.data.ensembles = cell(size(app.data,1),max_planes);

if strcmpi(app.RunallDropDown.Value, 'stats')
    fprintf('Running all stats, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        params.ddata = ddata;
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            if isempty(app.data(n_dset,:).stats{n_pl})
                params.cdata = app.data(n_dset,:).cdata{n_pl};
                app.data(n_dset,:).stats{n_pl} = f_dv_compute_stats_core(app, params);
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'data_dim_pca')
    fprintf('Running all dim red pca, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            if isempty(app.data(n_dset,:).data_dim_pca{n_pl})
                params.cdata = app.data(n_dset,:).cdata{n_pl};
                app.data(n_dset,:).data_dim_pca{n_pl} = f_dv_estimate_dim_pca_core(params);
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'data_dim_cv')
    fprintf('Running all dim red CV, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        params.ddata = ddata;
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            if isempty(app.data(n_dset,:).data_dim_cv{n_pl})
                params.cdata = app.data(n_dset,:).cdata{n_pl};
                if isempty(app.data(n_dset,:).data_dim_pca{n_pl})
                    app.data(n_dset,:).data_dim_pca{n_pl} = f_dv_estimate_dim_pca_core(params);
                    params.data_dim_pca = app.data(n_dset,:).data_dim_pca{n_pl};
                else
                    params.data_dim_pca = app.data(n_dset,:).data_dim_pca{n_pl};
                end
                app.data(n_dset,:).data_dim_cv{n_pl} = f_dv_estimate_dim_cv_core(params);
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'ensembles')
    fprintf('Running all ensemble extract, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        params.ddata = ddata;
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            if isempty(app.data(n_dset,:).ensembles{n_pl})
                params.cdata = app.data(n_dset,:).cdata{n_pl};
                if isempty(app.data(n_dset,:).data_dim_pca{n_pl})
                    app.data(n_dset,:).data_dim_pca{n_pl} = f_dv_estimate_dim_pca_core(params);
                    params.data_dim_pca = app.data(n_dset,:).data_dim_pca{n_pl};
                else
                    params.data_dim_pca = app.data(n_dset,:).data_dim_pca{n_pl};
                end
                app.data(n_dset,:).ensembles{n_pl} = f_dv_ensemble_extract_core(app, params);
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'ensemble_stats')
    fprintf('Running all ensemble stats, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        params.ddata = ddata;
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            if isempty(app.data(n_dset,:).ensemble_stats{n_pl})
                params.cdata = app.data(n_dset,:).cdata{n_pl};
                if isempty(app.data(n_dset,:).data_dim_pca{n_pl})
                    app.data(n_dset,:).data_dim_pca{n_pl} = f_dv_estimate_dim_pca_core(params);
                    params.data_dim_pca = app.data(n_dset,:).data_dim_pca{n_pl};
                else
                    params.data_dim_pca = app.data(n_dset,:).data_dim_pca{n_pl};
                end
                if isempty(app.data(n_dset,:).ensembles{n_pl})
                    app.data(n_dset,:).ensembles{n_pl} = f_dv_ensemble_extract_core(app, params);
                    params.ensembles = app.data(n_dset,:).ensembles{n_pl};
                else
                    params.ensembles = app.data(n_dset,:).ensembles{n_pl};
                end
                app.data(n_dset,:).ensemble_stats{n_pl} = f_dv_ensamble_stats_core(app, params);
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'ensemble_tuning')
    fprintf('Running all ensemble tuning, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        params.ddata = ddata;
        for n_pl = 1:ddata.num_planes
            params.n_pl = n_pl;
            if isempty(app.data(n_dset,:).ensemble_tuning{n_pl})
                params.cdata = app.data(n_dset,:).cdata{n_pl};
                if isempty(app.data(n_dset,:).data_dim_pca{n_pl})
                    app.data(n_dset,:).data_dim_pca{n_pl} = f_dv_estimate_dim_pca_core(params);
                    params.data_dim_pca = app.data(n_dset,:).data_dim_pca{n_pl};
                else
                    params.data_dim_pca = app.data(n_dset,:).data_dim_pca{n_pl};
                end
                if isempty(app.data(n_dset,:).ensembles{n_pl})
                    app.data(n_dset,:).ensembles{n_pl} = f_dv_ensemble_extract_core(app, params);
                    params.ensembles = app.data(n_dset,:).ensembles{n_pl};
                else
                    params.ensembles = app.data(n_dset,:).ensembles{n_pl};
                end
%                 if isempty(app.data(n_dset,:).ensemble_stats{n_pl})
%                     app.data(n_dset,:).ensemble_stats{n_pl} = f_dv_ensamble_stats_core(app, params);
%                     params.ensemble_stats = app.data(n_dset,:).ensemble_stats{n_pl};
%                 else
%                     params.ensemble_stats = app.data(n_dset,:).ensemble_stats{n_pl};
%                 end
                app.data(n_dset,:).ensemble_tuning{n_pl} = f_dv_ensemble_tuning(app, params);
            end
        end
    end
    fprintf('\n');
end
disp('Done')
end