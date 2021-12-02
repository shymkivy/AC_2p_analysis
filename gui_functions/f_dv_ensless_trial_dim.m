function f_dv_ensless_trial_dim(app)

% this will compute dimensionality of trials
params2 = f_dv_ensemble_params(app);
est_params_pca = params2.est_params_pca;

est_params_pca.dim_est_num_reps = 100;
%%

firing_rate = app.cdata.S;
trial_types = app.ddata.trial_types{1};
stim_times = app.ddata.stim_frame_index{1};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
mmn_freq = app.ddata.MMN_freq{1};
stats1 = app.ddata.stats{1};
ens_list = app.ddata.ensembles{1}.ens_out.cells.ens_list(app.ddata.ensemble_stats{1}.accepted_ensembles);

ens_cells = false(stats1.num_cells, numel(ens_list));
for n_ens = 1:numel(ens_list)
    ens_cells(ens_list{n_ens}, n_ens) = 1;
end

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trig_window);
[trial_data_sort_wctx, trial_types_wctx, trial_types_idx_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

%%

plot_individual_trials = app.PlotindivtrialsCheckBox.Value;
shuffle_trial_order = app.ShufletrorderinsteadCheckBox.Value;



tn_all = f_dv_get_trial_number(app);

num_tn = numel(tn_all);

dim_est_tr = zeros(num_tn,1);
num_cells_all = zeros(num_tn,1);
num_trials_all = zeros(num_tn,1);

hc_params.plot_dist_mat = 0;
hc_params.plot_clusters = 0;



for n_tt = 1:num_tn
    tn = tn_all(n_tt);
    tt = app.ops.context_types_all(tn);
    tr_idx = tt == trial_types_wctx;
    
    resp_cells = stats1.cell_is_resp(:,tn);
    
    resp_ens = logical(sum(app.ddata.ensemble_tuning{1}.cell_is_resp(:,tn),2));
    resp_ens_cells = ens_cells(:,resp_ens); 
    
    if app.PlotrespensCheckBox.Value
        resp_full = [resp_cells,resp_ens_cells];
    else
        resp_full = [resp_cells,ens_cells];
    end
    
    if strcmpi(app.SelectcellsDropDown.Value, 'Trial responsive')
        cell_idx = resp_cells;
    elseif strcmpi(app.SelectcellsDropDown.Value, 'Ensemble')
        cell_idx = logical(sum(resp_ens_cells,2));
    elseif strcmpi(app.SelectcellsDropDown.Value, 'Both')
        cell_idx = logical(sum([resp_cells, resp_ens_cells],2));
    elseif strcmpi(app.SelectcellsDropDown.Value, 'All cells')
        cell_idx = true(stats1.num_cells,1);
    end
    
    num_tr = sum(tr_idx);
    num_cells = sum(cell_idx);
    
    num_cells_all(n_tt) = num_cells;
    num_trials_all(n_tt) = num_tr;
    
    if num_cells
    
        data2 = trial_data_sort_wctx(cell_idx,:,tr_idx);
        data2_2d = reshape(data2,num_cells,[]);
        
        resp_full2 = resp_full(cell_idx,:);
        
        if num_cells > 1

            if shuffle_trial_order
                data_dim_est = f_ensemble_comp_data_dim_trials(data2, est_params_pca);
            else
                data_dim_est = f_ensemble_comp_data_dim2(data2_2d, est_params_pca);
            end

            hc_out_cell = f_hcluster_wrap(data2_2d, hc_params);

            data2_2d = data2_2d(hc_out_cell.dend_order,:);
            data2 = data2(hc_out_cell.dend_order,:,:);
            
            resp_full2 = resp_full2(hc_out_cell.dend_order,:);
            
            dim_est_tr(n_tt) = data_dim_est.dimensionality_corr;
        else
            dim_est_tr(n_tt) = 0;

            hc_out_cell.dend_order = 1;
        end
        
        hc_params_tr = hc_params;
        hc_params_tr.num_clust = app.NumclustEditField.Value;
        data2_2d_tr = reshape(data2,[], num_tr);
        hc_out_tr = f_hcluster_wrap(data2_2d_tr', hc_params_tr);
        
        clust_raster = false(num_tr, hc_out_tr.num_clust);
        for n_clust = 1:hc_out_tr.num_clust
            clust_raster(hc_out_tr.clust_ident == n_clust,n_clust) = 1;
        end
        clust_raster_sort = clust_raster(hc_out_tr.dend_order,:);
        
        data2_2d_tr_sort = reshape(data2(:,:,hc_out_tr.dend_order),num_cells,[]);
        dist1 = f_pdist_YS(data2_2d_tr', 'cosine');

        if plot_individual_trials
            figure; 
            subplot(9,11,[1:3 11*1+(1:3) 11*2+(1:3)]);
            imagesc(data2_2d);
            title(sprintf('%s, trial=%d, corr=%.2f', app.ddata.experiment{1}, tt, dim_est_tr(n_tt)), 'Interpreter', 'none');
            ylabel('Sorted cells')
            xlabel('frames');
            
            if num_cells > 1
                s1 = subplot(9,11,[4:6 11*1+(4:6) 11*2+(4:6)]);
                imagesc(1 - hc_out_cell.dist)
                sp = gca;
                sp.YDir = 'reverse';
                title('Cell - cell similarity')
                xlabel('cells');
                s1.YTick = [];
            end
            s1 = subplot(9,11,[7:9 11*1+(7:9) 11*2+(7:9)]);
            imagesc(data2_2d_tr_sort);
            title('Sorted trials');
            xlabel('frames');
            s1.YTick = [];
            
            s1 = subplot(9,11,[10:11 11*1+(10:11) 11*2+(10:11)]);
            imagesc(resp_full2)
            title('ensemble identity');
            xlabel('ensemble');
            s1.YTick = [];
            
            subplot(9,11, [11*4+(1:3) 11*5+(1:3) 11*6+(1:3)]);
            imagesc(1 - dist1)
            sp = gca;
            sp.YDir = 'reverse';
            title('Trial - trial similarity unsorted');
            xlabel('trials');
            
            subplot(9,11,11*8+(1:3));
            imagesc(clust_raster');
            title('trial cluster identity')
            xlabel('trials');
            
            subplot(9,11,[11*4+(7:9) 11*5+(7:9) 11*6+(7:9)]);
            imagesc(1 - hc_out_tr.dist)
            sp = gca;
            sp.YDir = 'reverse';
            title('Trial - trial similarity sorted');
            xlabel('trials');
            
            subplot(9,11,11*8+(7:9));
            imagesc(clust_raster_sort');
            title('sorted trial cluster identity')
            xlabel('trials');
            
        end
        
        if app.plotclusttriggeredstuffCheckBox.Value

            if strcmpi(app.clusttriggeredwhatDropDown.Value, 'stim')                
                num_frames = numel(app.ddata.proc_data{1}.stim_times_trace{1});
                num_freqs = app.ddata.proc_data{1}.stim_params.num_freqs;
                num_stim = numel(stim_times);
                stim_dur_frames = round(app.ddata.proc_data{1}.stim_params.stim_duration*1000/app.ddata.proc_data{1}.frame_data.volume_period);
                mmn_freq = app.ddata.MMN_freq{1};
                
                stim_trace = zeros(num_freqs, num_frames);
                
                for n_stim = 1:num_stim
                    stim_time = stim_times(n_stim);
                    stim_type = trial_types(n_stim);

                    if stim_type == 170
                        stim_type2 = mmn_freq(2);
                    elseif stim_type == 270
                        stim_type2 = mmn_freq(1);
                    elseif and(stim_type>100, stim_type<170)
                        stim_type2 = mmn_freq(1);
                    elseif and(stim_type>200, stim_type<270)
                        stim_type2 = mmn_freq(2);
                    else
                        stim_type2 = stim_type;
                    end
                    
                    stim_trace(stim_type2, stim_time:(stim_time+stim_dur_frames)) = 1;
                      
                end
                trial_list = trial_types_idx_wctx(tr_idx);
                trial_list_clust = trial_list(hc_out_tr.clust_ident == app.ClustertouseEditField.Value);
                
                stim_times_trial = stim_times(trial_list);
                stim_times_clust = stim_times(trial_list_clust);
                stim_trig_window = [27 27];
                num_bins = sum(stim_trig_window);
                
                stim_trial_data_sort = f_get_stim_trig_resp(stim_trace, stim_times_trial, stim_trig_window);
                trial_trig_ave = mean(stim_trial_data_sort,3);
                
                stim_clust_data_sort = f_get_stim_trig_resp(stim_trace, stim_times_clust, stim_trig_window);
                clust_trig_ave = mean(stim_clust_data_sort,3);
                
                figure; hold on; axis tight;
                subplot(1,10,1:8);
                imagesc(reshape(stim_clust_data_sort, num_freqs, []));
                for n_stim = 1:numel(stim_times_clust)
                    line([n_stim n_stim]*num_bins, [0.5 num_freqs+0.5], 'color','r');
                end
                title(sprintf('%s, trial=%d, corr=%.2f, clust=%d', app.ddata.experiment{1}, tt, dim_est_tr(n_tt),app.ClustertouseEditField.Value), 'Interpreter', 'none');
                subplot(1,10,9);
                imagesc(clust_trig_ave);
                title('clust trig ave');
                subplot(1,10,10);
                imagesc(trial_trig_ave)
                title('trial trig ave');
                
            elseif strcmpi(app.clusttriggeredwhatDropDown.Value, 'ensembles')
            elseif strcmpi(app.clusttriggeredwhatDropDown.Value, 'cells')
            end
            
        end
        
    else
        dim_est_tr(n_tt) = 0;
    end
end


if app.PlotcorrtotalsCheckBox.Value
    c1 = categorical(app.ops.context_types_labels(tn_all));
    c2 = reordercats(c1,app.ops.context_types_labels(tn_all));

    figure; 
    subplot(2,1,1);
    bar(c2,dim_est_tr)
    title(sprintf('%s', app.ddata.experiment{1}), 'Interpreter', 'none');
    subplot(2,1,2);
    bar(c2,num_cells_all)
    title('Num cells');
end


end