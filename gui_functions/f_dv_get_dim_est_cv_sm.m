function array_out = f_dv_get_dim_est_cv_sm(app)

cent1 = app.dimestcv_sm_centerEditField.Value;
range1 = app.dimestcv_sm_rangeEditField.Value;
count1 = app.dimestcv_sm_countEditField.Value;

shift = (min(cent1-range1,0));

array_out = round(linspace(cent1-range1-shift, cent1+range1-shift, count1));


end