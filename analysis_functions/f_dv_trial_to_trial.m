function f_dv_trial_to_trial(app)

% app stuff
ops = app.ops;
params = f_dv_gather_params(app);

tn_all = f_dv_get_trial_number(params);
[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);

% dset_list = 1:3;
% dset_list = 1:7;

do_mean = 1;
sort_plot = 'clust';   % none, clust, hist

num_dsets = size(data,1);

if do_mean
    title_tag2 = sprintf('%s; mean; sort %s', title_tag, sort_plot);
    clim2 = [-0.4, 1];
else
    title_tag2 = sprintf('%s; full; sort %s', title_tag, sort_plot);
    clim2 = [-0.2, 0.8];
end

%tn_all = 2:9;
num_tn = numel(tn_all);

% select_resp_cells = app.selectrespcellsCheckBox.Value;
% resort_by_ens = app.resortbyensCheckBox.Value;
% sort_trials = app.sorttrialsCheckBox.Value;
% sort_with_full_firing_rate = app.sortwithfullfrCheckBox.Value;

corr_vals = nan(num_dsets, num_tn);
isi_vals = nan(num_dsets,1);
use_dset = false(num_dsets, 1);
for n_dset = 1:num_dsets
    %dset_idx = dset_list(n_dset);
    ddata = data(n_dset,:);
    isi_vals(n_dset) = round(median(diff(ddata.proc_data{1}.stim_times_volt{1}))/1000 - 0.5,1);
    
    if strcmpi(ddata.paradigm, 'cont')
        use_dset(n_dset) = 1;
    end
    
    if use_dset(n_dset)

        cdata = ddata.cdata{params.planes};
        stats1 = ddata.stats{params.planes};

        %num_cells = sum([cdata.num_cells]);
        firing_rate = cat(1,cdata.S_sm);
        trial_types = ddata.trial_types{1};
        trial_types_hist = [0; trial_types(1:end-1)];
        
        stim_frame_index = ddata.stim_frame_index{1};
        
        [~, trial_frames] = f_dv_compute_window_t(params.trial_window, mean(cat(1,cdata.volume_period)));
        
        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_frames);

        [selected_cells, ~, ~, ~, resp_cells] = f_dv_get_resp_vals_cells(stats1, tn_all, params);

        for n_tn = 1:num_tn
            tn1 = tn_all(n_tn);
            
            tr_idx = logical(sum(trial_types == ops.context_types_all(tn1)',2));
            tr_data = trial_data_sort(:,:,tr_idx);
            tr_hist = trial_types_hist(tr_idx);
            
            num_trials = sum(tr_idx);

            selected_cells2 = selected_cells(:,n_tn);
            resp_cells2 = resp_cells(:,n_tn);
            
            if sum(resp_cells2)>5
                
                tr_data2 = tr_data(selected_cells2,:,:);
                %firing_rate2 = firing_rate(resp_cells,:);
                %tr_data3 = squeeze(mean(tr_data2,2));

                hc_params.plot_dist_mat = 0;
                hc_params.plot_clusters = 0;
                hc_params.num_clust = 1;
        
                if do_mean
                    tr_data_2d_tr = squeeze(mean(tr_data2,2));
                else
                    tr_data_2d_tr = reshape(tr_data2, [], num_trials);
                end

                hclust_out_trial = f_hcluster_wrap(tr_data_2d_tr', hc_params);
                
                SI = 1 - hclust_out_trial.dist_no;

                if strcmpi(sort_plot, 'clust')
                    sort_idx = hclust_out_trial.dend_order;
                    
                    SI = SI(sort_idx,:);
                    SI = SI(:, sort_idx);

                elseif strcmpi(sort_plot, 'hist')
                    [~, sort_idx] = sort(tr_hist);

                    SI = SI(sort_idx,:);
                    SI = SI(:, sort_idx);
                end

                figure;
                imagesc(SI);
                axis equal tight
                xlabel('trials');
                ylabel('trials');
                title(sprintf('trial %d; isi %.1fsec; %s', tn1, isi_vals(n_dset), title_tag2))
                clim(clim2);

            end
       
        end
    end
end



end
