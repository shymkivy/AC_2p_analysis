function f_dv_load_reg_data(app)

disp('Loading reg data...');
app.reg_data = load(app.regdatapathEditField.Value);
disp('Done');

end