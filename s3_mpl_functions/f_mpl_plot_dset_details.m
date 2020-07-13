function f_mpl_plot_dset_details(data, ops)
%% for each dset plot things

sp_h = {};
ylim1 = [0, 0.2];

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    for n_dset = 1:cdata.num_dsets
        
        if ops.use_zscores
            trial_ave = cdata.trial_ave_z{n_dset};
        else
            trial_ave = cdata.trial_ave{n_dset};
        end
        
        tuning_all = cdata.tuning_all{n_dset};
        ctx_mmn = cdata.ctx_mmn{n_dset};
        resp_temp = logical(tuning_all.peak_tuned_trials_combined_ctx);
        
        %% Plot all ctx
        if ops.plot.ctx_full_each_dset
            f_mpl_plot_ctx2(trial_ave, resp_temp, ctx_mmn, cdata.trial_window_t{n_dset}, ops);
            suptitle(sprintf('%s, dset %d', cond_name, n_dset));
        end
        
        %% Plot dd plot
        if ops.plot.dd_each_dset
            [sp_h1, ylim2] = f_mpl_plot_dd(trial_ave(:,:,ctx_mmn), resp_temp,cdata.trial_window_t{n_dset}, ops);
            suptitle(sprintf('%s, dset %d', cond_name, n_dset));
            sp_h = cat(1,sp_h, sp_h1);
            ylim1 = [min([ylim1(1), ylim2(1)]) max([ylim1(2), ylim2(2)])];
        end
        
        if ops.spat_tuning_plot
            %f_mpl_plot_spatial_tunning(cdata, ops, cond_name, n_dset);
        end
        
    end
end

% make all limits equal or dd plot
for n_pl = 1:numel(sp_h)
    ylim(sp_h{n_pl}, ylim1);
end






end