function f_mpl_population_analysis_raster(data, ops) 
%% input parameters for cross validation estimation of smooth window and number of correlated components / ensembles

estimate_params = 0;    % do estimation?
est_params.ensamble_method = 'nmf';              % options: svd, nmf, ica                % SVD is most optimal for encoding, NMF rotates components into something that is real and interpretable
est_params.normalize = 'norm_mean_std'; % 'norm_mean_std', 'norm_mean' 'none'   % either way, need to normalize the power of signal in each cell, otherwise dimred will pull out individual cells
est_params.smooth_SD = 120;       % range of values to estimate across    % larger window will capture 'sequences' of ensembles, if window is smaller than optimal, you will end up splitting those into more components
est_params.num_comp = 10:2:20;               % range of values to estimate across    
est_params.shuffle_data_chunks = 1;   % 1 or 0, keeping cell correlations   % if the sequence of trial presentation contains information, you will need to shuffle. Also need to do in chunks because adjacent time bins are slightly correlated
est_params.reps = 2;                   % how many repeats per param 

%%
est_params.n_rep = 1:est_params.reps;
est_params_list = f_build_param_list(est_params, {'smooth_SD', 'num_comp', 'n_rep'});

%% input paramseters for ensemble analysis
% NMF ensemble detection is best with thresh extraction
% for NMF best to use norm_rms(keep values positive), otherwise can also use norm_mean_std
% NMF 14 comp
% SVD 11-14 comp?
ens_params.ensamble_method = 'nmf'; % options: svd, nmf, ica     % here NMF is
ens_params.num_comp = 15;
ens_params.smooth_SD = 120; % 110 is better?
ens_params.normalize = 'norm_mean_std'; % 'norm_mean_std', 'norm_mean' 'none'
ens_params.ensamble_extraction = 'thresh'; %  'thresh'(for nmf) 'clust'(for svd)
ens_params.ensamble_extraction_thresh = 'signal_z'; % 'shuff' 'signal_z' 'signal_clust_thresh'
ens_params.signal_z_thresh = 2.5;
ens_params.shuff_thresh_percent = 95;
ens_params.plot_stuff = 0;

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
                
                firing_rate_norm = f_normalize(firing_rate, est_params.normalize);
                
                params1 = est_params_list(n_par);
                params1.vol_period = vol_period;
                accuracy = f_ens_estimate_corr_dim_cv(firing_rate_norm, params1);
                
                temp_fields = fields(accuracy);
                for n_fl = 1:numel(temp_fields)
                    est_params_list(n_par).(temp_fields{n_fl}) = accuracy.(temp_fields{n_fl});
                end
            end
            fprintf('\nDone\n');
            [~, min_ind] = min([est_params_list.test_err]);
            fprintf('From provided range, optimal smooth_SD = %d; Number of CV %s num_comp = %d\n', est_params_list(min_ind).smooth_SD, est_params.method, est_params_list(min_ind).num_comp);
            
            f_plot_cv_error_3D(est_params_list, 'smooth_SD', 'num_comp', 'test_err');
            ax1 = gca;
            ax1.Title.String = sprintf('%s, dset%d; %s L2 error from raw, (%s)',cond_name,n_dset,est_params.method, ax1.Title.String);          
     end

        %% Smooth data
        firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/vol_period);
        
        %% extract ensambles
        ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm, ens_params);
        
        %% analyze ensembles
        f_plot_raster_mean(firing_rate_sm(ens_out.ord_cell,:), 1);
        title('raster cell sorted');
        
        for n_comp = 1:numel(ens_out.cells.ens_list)
            cells1 = ens_out.cells.ens_list{n_comp};
            trials1 = ens_out.trials.ens_list{n_comp};
            scores1 = ens_out.scores(ens_out.cells.scores_alignment(n_comp),:);
            
            f_plot_ensamble_deets(firing_rate_sm, cells1, trials1, scores1);
            title([ens_params.ensamble_method ' ensamble ' num2str(n_comp)]);
        end
        
        
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


