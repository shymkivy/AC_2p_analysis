function f_mpl_plot_cond_details(data, ops)

%%
f_mpl_plot_dd_cond(data, ops);

%%
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    
    if ops.use_zscores
        trial_ave = cat(1,cdata.trial_ave_z{:});
    else
        trial_ave = cat(1,cdata.trial_ave{:});
    end
    
    ctx_mmn_full = cell(numel(cdata.area),1);
    resp_cells = cell(numel(cdata.area),1);
    for n_dset = 1:numel(cdata.area)
        ctx_mmn_full{n_dset} = repmat(cdata.ctx_mmn{n_dset}, cdata.num_cells(n_dset),1);
        resp_cells{n_dset} = cdata.peak_tuned_trials_combined_ctx{n_dset};
    end
    ctx_mmn_full = cat(1,ctx_mmn_full{:});
    resp_cells = cat(1,resp_cells{:});
    
    %%
    if sum(resp_cells(:))
        f_mpl_plot_ctx3(trial_ave, resp_cells, cdata.trial_window{1}.trial_window_t, ops);
        suptitle(sprintf('%s', cond_name));
    end
    %f_mpl_plot_ctx2(trial_ave, resp_cells,ctx_mmn_full, cdata.trial_window{1}.trial_window_t, ops);
    
end

% if ops.ctx_plots
%     f_mpl_plot_ctx_cond(data, ops);
% end

%%
if ops.tuning_plots
    f_mpl_plot_tuning_cond(data, ops);
end



end