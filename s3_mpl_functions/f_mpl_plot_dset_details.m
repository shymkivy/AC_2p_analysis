function f_mpl_plot_dset_details(data, ops)
%% for each dset plot things

sp_h = {};
ylim1 = [0, 0.2];

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    
    for n_dset = 1:numel(cdata.area)
        if 1 %cdata.num_cells(n_dset) == max(cdata.num_cells)
        
            if ops.use_zscores
                trial_ave = cdata.trial_ave_z{n_dset};
            else
                trial_ave = cdata.trial_ave{n_dset};
            end
            trial_window_t = cdata.trial_window{n_dset}.trial_window_t;

            ctx_mmn = cdata.ctx_mmn{n_dset};
            resp_cells_ctx = logical(cdata.peak_tuned_trials_combined_ctx{n_dset});
            
            dset_params.ctx_mmn = ctx_mmn;
            dset_params.cond_name = cond_name;
            dset_params.n_dset = n_dset;
            
            %% Plot all ctx
            if ops.plot.ctx_full_each_dset
                f_mpl_plot_ctx3(trial_ave, resp_cells_ctx, trial_window_t, ops);
                %f_mpl_plot_ctx4(trial_ave, resp_cells_ctx, trial_window_t, ops);
                suptitle(sprintf('%s, dset %d', cond_name, n_dset));
            end

            %% Plot dd plot
            if ops.plot.dd_each_dset
                [sp_h1, ylim2] = f_mpl_plot_dd(trial_ave(:,:,ctx_mmn), resp_cells_ctx,trial_window_t, ops);
                suptitle(sprintf('%s, dset %d', cond_name, n_dset));
                sp_h = cat(1,sp_h, sp_h1);
                ylim1 = [min([ylim1(1), ylim2(1)]) max([ylim1(2), ylim2(2)])];
            end

            if ops.plot.spatial_tuning_dset
                f_mpl_plot_spatial_tunning(cdata, ops, cond_name, n_dset);
            end
            
            if ops.plot.tuning_dset
                resp_cells_on = logical(cdata.peak_tuned_trials_onset{n_dset});
                resp_cells_off = logical(cdata.peak_tuned_trials_offset{n_dset});
                f_mpl_plot_tuning_dset(resp_cells_on, resp_cells_off, dset_params, ops);
            end
        end
            
    end
end

% make all limits equal or dd plot
for n_pl = 1:numel(sp_h)
    ylim(sp_h{n_pl}, ylim1);
end






end