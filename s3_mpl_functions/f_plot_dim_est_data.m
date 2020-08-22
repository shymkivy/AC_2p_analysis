function f_plot_dim_est_data(dim_est_st, ops)

max_cells = 30;

f_plot_dset_param_v_ncell(dim_est_st, 'dimensionality_corr', max_cells, ops);

f_plot_dset_param_v_ncell(dim_est_st, 'dimensionality_total', max_cells, ops);

f_plot_dset_param_v_ncell(dim_est_st, 'dimensionality_total_norm', max_cells, ops);

f_plot_dset_param_v_ncell(dim_est_st, 'dimensionality_total_norm_shuff', max_cells, ops);

f_plot_dset_param_v_ncell(dim_est_st, 'dimensionality_first_comp_size', max_cells, ops);
ylabel('percent variance')

end