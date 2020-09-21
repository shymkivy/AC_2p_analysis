function f_mpl_population_analysis_raster(data, ops) 
%%
estimate_params = 0;
est_params.method = 'svd'; % options: svd, nmf, ica
est_params.normalize = 'norm_rms'; % 'norm_mean_std', 'norm_mean' 'none'
est_params.smooth_SD = 100:100:500;   % range of values to estimate across
est_params.num_comp = 14;       % range of values to estimate across
est_params.randomize_trials = 0;
est_params.n_rep = 1;
% 100 best with none or norm_mean
% 0 > 35 > 50 with norm_full
est_params_list = f_build_param_list(est_params, {'smooth_SD', 'num_comp', 'n_rep'});

%%
% NMF ensemble detection is best
% for NMF best to use norm_rms(keep values positive), otherwise can also use norm_mean_std
% NMF 14 comp
% SVD 11-14 comp?
ens_params.method = 'NMF'; % options: svd, nmf, ica
ens_params.num_comp = 14;
ens_params.smooth_SD = 100; % 110 is better?
ens_params.normalize = 'norm_rms'; % 'norm_mean_std', 'norm_mean' 'none'
ens_params.ensamble_extraction = 'thresh'; % 'clust'(for svd) 'thresh'(for nmf)
ens_params.ensamble_extraction_thresh = 'shuff'; % 'signal_z' 'shuff' 'signal_clust_thresh'
ens_params.plot_stuff = 1;

%%

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
        if estimate_params
            fprintf('Estimating params n/%d reps: ',numel(est_params_list));
            %dim_corr = zeros(numel(estimate_smooth_list),1);
            for n_par = 1:numel(est_params_list)
                fprintf('%d..',n_par);
                
                firing_rate_norm = f_normalize(firing_rate, ens_params.normalize);
                
                params1 = est_params_list(n_par);
                params1.vol_period = vol_period;
                accuracy = f_ens_estimate_corr_dim_cv(firing_rate_norm, params1);
                
                est_params_list(n_par).train_err = accuracy.train_err;
                est_params_list(n_par).train_err_sm = accuracy.train_err_sm;
                est_params_list(n_par).test_err = accuracy.test_err;
                est_params_list(n_par).test_err_sm = accuracy.test_err_sm;
            end
            fprintf('\nDone\n');
            [~, min_ind] = min([est_params_list.test_err_sm]);
            fprintf('From provided range, optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', est_params_list(min_ind).smooth_SD, est_params.method, est_params_list(min_ind).num_comp);
            
            f_plot_cv_error_3D(est_params_list, 'smooth_SD', 'num_comp', 'test_err');
            ax1 = gca;
            ax1.Title.String = sprintf('%s, dset%d; %s error from raw, %s',cond_name,n_dset,est_params.method, ax1.Title.String);
            
            f_plot_cv_error_3D(est_params_list, 'smooth_SD', 'num_comp', 'test_err_sm');
            ax1 = gca;
            ax1.Title.String = sprintf('%s, dset%d; %s error from smooth, %s',cond_name,n_dset,est_params.method, ax1.Title.String);
        end

        %% Smooth data
        firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/vol_period);
        
        %%
        ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm, ens_params);
        
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


