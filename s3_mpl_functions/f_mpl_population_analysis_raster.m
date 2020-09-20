function f_mpl_population_analysis_raster(data, ops) 
estimate_smooth = 0;
estimate_smooth_list = 0:50:500;
% 100 best with none or norm_mean
% 0 > 35 > 50 with norm_full
smooth_sd = 100;

norm_method = 'norm_full'; % 'norm_full', 'norm_mean' 'none'

disp('Ensemble analysis...');
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    for n_dset = 1:cdata.num_dsets
        
        firing_rate = cat(1,cdata.firing_rate{n_dset,:});
        mmn_phase_mpl = cdata.proc_data{n_dset}.mmn_phase_mpl{1};
        firing_rate_cont = firing_rate(:,mmn_phase_mpl == 1);
        vol_period = cdata.proc_data{n_dset}.frame_data.volume_period;
        
        %firing_rate = firing_rate_smooth(:,sum(firing_rate_smooth) >0);
        firing_rate = firing_rate_cont;
        
        %% remove inactive cells
        
        active_cells = sum(firing_rate,2) > 0;
        firing_rate(~active_cells,:) = [];
        
        num_cells = size(firing_rate,1);
        
        firing_rate = firing_rate(randperm(num_cells),:);
        %% estimate best smoothing window
        %% estimate smooth
        if estimate_smooth
            fprintf('Estimating smooth SD n/%d reps: ',numel(estimate_smooth_list));
            dim_corr = zeros(numel(estimate_smooth_list),1);
            for n_sm = 1:numel(estimate_smooth_list)
                fprintf('%d..',n_sm);
                if estimate_smooth_list(n_sm)>0
                    raster_sm = f_smooth_gauss(firing_rate, estimate_smooth_list(n_sm)/vol_period);
                else
                    raster_sm = firing_rate;
                end
                firing_rate_sm = f_normalize(raster_sm, norm_method);
                
                dim_corr(n_sm) = f_ens_estimate_dim(firing_rate_sm);
            end
            figure; plot(estimate_smooth_list, dim_corr);
            xlabel('Smooth SD'); ylabel('Dimensionality of corr');
            title(sprintf('%s smooth', norm_method));
            fprintf(' Done\n');
        end

        %% Smooth data
        if smooth_sd> 0
            firing_rate_sm = f_smooth_gauss(firing_rate, smooth_sd/vol_period);
        else
            firing_rate_sm = firing_rate;
        end
        
        %%
        ens_params.plot_stuff = 1;
        ens_params.sort_cells = 1;
        ens_params.normalize = norm_method; % 'norm_full', 'norm_mean' 'none'
        ens_params.ensamble_extraction = 'thresh'; % clust 'thresh'
        ens_params.ensamble_extraction_thresh = 'shuff'; % 'signal_z' 'shuff' 'signal_clust_thresh'
        ens_out = f_ensemble_analysis_YS_raster(firing_rate, ens_params);
        
        %% analyze ensembles
        1
        
        %% ensemble analysis with Luis method
        
%         binary_firing_rate_smooth = double(firing_rate>2*std(firing_rate,[],2));
%         bin_no_act = sum(binary_firing_rate_smooth,2) == 0;
%         binary_firing_rate_smooth(bin_no_act,:) = [];
%         
%         % first make minary data
%         addpath('C:\Users\rylab_dataPC\Desktop\Yuriy\SVD_ensemble_Luis\SVDEnsemble-master');
%         params.pks = 4;
%         params.ticut = 0.2;
%         params.jcut = 0.06;
%         [core_svd,state_pks_full,param] = findSVDensemble(binary_firing_rate_smooth,[],params);
%         
%         figure; hold on;
%         temp_trace = core_svd{3};
%         for n_cell = 1
%             plot(sum(firing_rate(temp_trace,:)))
%         end
%         figure; hold on;
%         plot(firing_rate_smooth(10,:)*5)
%         plot(a(10,:))
        
    end
end


end


