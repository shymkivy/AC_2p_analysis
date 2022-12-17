function f_dv_plot_pop_mmn(app)

n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

ctx_plot_list = [18, 19, 20; ...
                 28, 29, 30]';

red_start = [201, 101];
             
             
num_dsets = numel(data.dset_name);
num_flips = 2;

params = f_dv_gather_params(app);

num_red = 20;

data_dset = cell(num_dsets,1);
vol_per_all = zeros(num_dsets,1);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    stats1 = ddata.stats{n_pl};
    params.n_dset = find(ddata.idx == app.data.idx);
    cdata = f_dv_compute_cdata(ddata, params);

    firing_rate = cdata.S_sm;
    trial_types = ddata.trial_types{1};
    stim_times = ddata.stim_frame_index{n_pl};
    mmn_freq = ddata.MMN_freq{1};
    
    data_temp = cell(2,1);
    
    for n_flip = 1:num_flips
        tn_all = ctx_plot_list(:,n_flip);
    
        sel_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all);
        sel_cells2 = logical(sum(sel_cells,2));
        
        mmn_start = find(trial_types == red_start(n_flip), 1);
        
        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times(mmn_start:(mmn_start+num_red-1)), trial_frames);
        
        trial_data_sort2 = trial_data_sort(sel_cells2,:,:);
        
        trial_data_sort3 = squeeze(mean(trial_data_sort2,1));
        
        data_temp{n_flip} = trial_data_sort2;
    end
    
    data_dset{n_dset} = cat(1, data_temp{:});
end

data_dset2 = cat(1, data_dset{:});

mean_data = squeeze(mean(data_dset2,1));

%ylim1 = [min(trial_data_sort3(:)), max(trial_data_sort3(:))];
figure;
plot(mean_data(:))

%%

data_dset = cell(num_dsets,1);
vol_per_all = zeros(num_dsets,1);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    stats1 = ddata.stats{n_pl};
    params.n_dset = find(ddata.idx == app.data.idx);
    cdata = f_dv_compute_cdata(ddata, params);

    firing_rate = cdata.S_sm;
    trial_types = ddata.trial_types{1};
    stim_times = ddata.stim_frame_index{n_pl};
    mmn_freq = ddata.MMN_freq{1};
    
    data_temp = cell(2,1);
    
    for n_flip = 1:num_flips
        tn_all = ctx_plot_list(:,n_flip);
    
        sel_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all);
        sel_cells2 = logical(sum(sel_cells,2));
        
        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);

        trial_ave = f_mpl_trial_average(trial_data_sort, trial_types, red_start(n_flip):(red_start(n_flip)+7), 'none');
        
        trial_data_sort2 = trial_ave(sel_cells2,:,:);
        
        data_temp{n_flip} = trial_data_sort2;
    end
    
    data_dset{n_dset} = cat(1, data_temp{:});
end

data_dset2 = cat(1, data_dset{:});


mean_data = squeeze(mean(data_dset2,1));

%ylim1 = [min(trial_data_sort3(:)), max(trial_data_sort3(:))];
figure;
plot(mean_data(:))
ylim([0 .1])


end