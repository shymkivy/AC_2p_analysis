function cdata = f_dv_get_cdata(app)

if strcmpi(app.SelectdatagroupButtonGroup.SelectedObject.Text, 'plane')
    n_pl = app.mplSpinner.Value;
    cdata = app.cdata{n_pl};
else
    cdata = cat(1,app.cdata{:});
end

end