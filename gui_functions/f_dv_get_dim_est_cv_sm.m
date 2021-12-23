function array_out = f_dv_get_dim_est_cv_sm(app)

cent1 = app.dimestcv_sm_centerEditField.Value;
range1 = app.dimestcv_sm_rangeEditField.Value;
count1 = app.dimestcv_sm_countEditField.Value;

array_out = round(linspace(max(cent1-range1, 0), cent1+range1, count1));


end