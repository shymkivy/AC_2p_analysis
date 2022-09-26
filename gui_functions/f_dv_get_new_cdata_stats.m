function [cdata, stats] = f_dv_get_new_cdata_stats(app, ddata, params)

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    cdata = f_dv_compute_cdata(ddata, params);
    stats = ddata.stats{app.mplSpinner.Value};
else
    cdata = cell(ddata.num_planes,1);
    for n_pl = 1:ddata.num_planes
        params.n_pl = n_pl;
        cdata{n_pl} = f_dv_compute_cdata(ddata, params);
    end
    cdata = cat(1,cdata{:});
    stats = cat(1,ddata.stats{:});
end

end