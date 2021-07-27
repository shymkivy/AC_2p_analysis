function f_dv_plot_ens_counts(app)


ens_resp = app.ddata.ensemble_tuning{1}.cell_is_resp;

figure;
bar(sum(ens_resp))



end