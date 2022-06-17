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
            if or(isempty(ddata.stats{n_pl}), app.OverwriteCheckBox.Value)
                params.cdata = f_dv_compute_cdata(app, params);
                if isfield(params.ddata.proc_data{1}, 'stim_params')
                    app.data(n_dset,:).stats{n_pl} = f_dv_compute_stats_core(app, params);
                else
                    fprintf('skipping %s; no stim_params\n', params.ddata.dset_name_full{1});
                end
            end
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'data_dim_pca')
    fprintf('Running all dim red pca, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        params.ddata = app.data(n_dset,:);
        if or(isempty(app.data(n_dset,:).data_dim_pca{1}), app.OverwriteCheckBox.Value)
            params.cdata = cat(1,params.ddata.cdata{:});
            app.data(n_dset,:).data_dim_pca{1} = f_dv_estimate_dim_pca_core(params);
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
        if or(isempty(ddata.data_dim_cv{1}), app.OverwriteCheckBox.Value)
            params.cdata = cat(1,params.ddata.cdata{:});
            if isempty(ddata.data_dim_pca{1})
                app.data(n_dset,:).data_dim_pca{1} = f_dv_estimate_dim_pca_core(params);
            end
            params.data_dim_pca = app.data(n_dset,:).data_dim_pca{1};
            app.data(n_dset,:).data_dim_cv{1} = f_dv_estimate_dim_cv_core(params);
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
        if or(isempty(app.data(n_dset,:).ensembles{1}), app.OverwriteCheckBox.Value)
            params.cdata = cat(1,params.ddata.cdata{:});
            if isempty(app.data(n_dset,:).data_dim_pca{1})
                app.data(n_dset,:).data_dim_pca{1} = f_dv_estimate_dim_pca_core(params);
            end
            params.data_dim_pca = app.data(n_dset,:).data_dim_pca{1};
            app.data(n_dset,:).ensembles{1} = f_dv_ensemble_extract_core(params);
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
        if or(isempty(app.data(n_dset,:).ensemble_stats{1}), app.OverwriteCheckBox.Value)
            params.cdata = cat(1,params.ddata.cdata{:});
            if isempty(app.data(n_dset,:).data_dim_pca{1})
                app.data(n_dset,:).data_dim_pca{1} = f_dv_estimate_dim_pca_core(params);
            end
            params.data_dim_pca = app.data(n_dset,:).data_dim_pca{1};
            if isempty(app.data(n_dset,:).ensembles{1})
                app.data(n_dset,:).ensembles{1} = f_dv_ensemble_extract_core(params);
            end
            params.ensembles = app.data(n_dset,:).ensembles{1};
            app.data(n_dset,:).ensemble_stats{1} = f_dv_ensemble_stats_core(params);
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'ensemble_tuning_stats')
    fprintf('Running all ensemble tuning, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        params.ddata = ddata;
        if or(isempty(app.data(n_dset,:).ensemble_tuning_stats{1}), app.OverwriteCheckBox.Value)
            params.cdata = cat(1,params.ddata.cdata{:});
            if isempty(app.data(n_dset,:).data_dim_pca{1})
                app.data(n_dset,:).data_dim_pca{1} = f_dv_estimate_dim_pca_core(params);
            end
            params.data_dim_pca = app.data(n_dset,:).data_dim_pca{1};
            if isempty(app.data(n_dset,:).ensembles{1})
                app.data(n_dset,:).ensembles{1} = f_dv_ensemble_extract_core(params);
            end
            if ~isempty(app.data(n_dset,:).ensemble_stats{1})
                params.ensemble_stats = app.data(n_dset,:).ensemble_stats{1};
            end
            params.ensembles = app.data(n_dset,:).ensembles{1};
            app.data(n_dset,:).ensemble_tuning_stats{1} = f_dv_ensemble_tuning_stats2(app, params);
        end
    end
    fprintf('\n');
elseif strcmpi(app.RunallDropDown.Value, 'ensless_dim_est')
    fprintf('Running all ensless dim est, Dset_/%d: ', num_data)
    for n_dset = 1:num_data
        fprintf('%d..', n_dset);
        params.n_dset = n_dset;
        ddata = app.data(n_dset,:);
        params.ddata = ddata;
        if or(isempty(app.data(n_dset,:).ensless_dim_est{1}), app.OverwriteCheckBox.Value)
            params.cdata = cat(1,params.ddata.cdata{:});
            app.data(n_dset,:).ensless_dim_est{1} = f_dv_ensless_dim_est(app, params);
        end
    end
    fprintf('\n');
end
disp('Done')
end