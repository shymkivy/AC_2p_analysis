function array_out = f_dv_get_dim_est_cv_nc(app)

cent1 = app.dimestcv_nc_centerEditField.Value;

if app.dimestcv_centeratdimpcaCheckBox.Value
    if app.DimpcaEditField.Value>0
        cent1 = app.DimpcaEditField.Value;
    end
end

range1 = app.dimestcv_nc_rangeEditField.Value;
count1 = app.dimestcv_nc_countEditField.Value;

shift = (min(cent1-range1,1)-1);

array_out = unique(round((linspace(cent1-range1+shift, cent1+range1+shift, count1))));

end