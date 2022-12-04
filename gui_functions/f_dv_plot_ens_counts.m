function f_dv_plot_ens_counts(app)

ens_resp = app.ddata.ensemble_tuning_stats{1}.resp_cells_peak;

figure;
bar(sum(ens_resp))

end