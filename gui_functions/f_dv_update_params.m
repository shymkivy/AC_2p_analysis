function f_dv_update_params(app)

if ~isempty(app.cp_app)
    if isvalid(app.cp_app)
        ddata = app.ddata;
        n_pl = app.mplSpinner.Value;
        heading1 = {'Dataset', ddata.experiment{1}; 'Plane', num2str(n_pl)};
        fieldnames1 = [];
        vals = [];
        if strcmpi(app.cp_app.SourceDropDown.Value, 'stats_params')
            stats1 = ddata.stats{n_pl};
            if ~isempty(stats1)
                params1 = stats1.stat_params;
                fieldnames1 = fieldnames(params1);
                vals = cell(numel(fieldnames1),1);
                for n_fl = 1:numel(fieldnames1)
                    vals{n_fl} = num2str(params1.(fieldnames1{n_fl}));
                end
            end
            
        elseif strcmpi(app.cp_app.SourceDropDown.Value, 'est_dim_pca')
            data_dim_pca = ddata.data_dim_pca{1};
            if ~isempty(data_dim_pca)
                params1 = data_dim_pca.params;
                if isfield(params1, 'params_in')
                    params1 = rmfield(params1, 'params_in');
                end
                fieldnames1 = fieldnames(params1);
                vals = cell(numel(fieldnames1),1);
                for n_fl = 1:numel(fieldnames1)
                    vals{n_fl} = num2str(params1.(fieldnames1{n_fl}));
                end
            end
        elseif strcmpi(app.cp_app.SourceDropDown.Value, 'est_dim_cv')
            data_dim_cv = ddata.data_dim_cv{1};
            if ~isempty(data_dim_cv)
            end
        elseif strcmpi(app.cp_app.SourceDropDown.Value, 'ens_params')
            data_ensembles = ddata.ensembles{1};
            if ~isempty(data_ensembles)
            else
            end
        end
        app.cp_app.UITable.Data = table([heading1; fieldnames1, vals]);
    end
end

end