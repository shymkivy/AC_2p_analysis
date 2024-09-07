function f_dv_ensless_single_trial_corr(app)

n_pl = app.mplSpinner.Value;

% dset_list = 1:3;
% dset_list = 1:7;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = size(data,1);

%tn_all = 2:9;
tn_all = f_dv_get_trial_number(app);
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

        if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
            cdata = ddata.cdata{n_pl};
            stats1 = ddata.stats{n_pl};
        else
            cdata = cat(1,ddata.cdata{:});
            stats1 = cat(1,ddata.stats{:});
        end

        %num_cells = sum([cdata.num_cells]);
        firing_rate = cat(1,cdata.S_sm);
        trial_types = ddata.trial_types{1};
        
        stim_frame_index = ddata.stim_frame_index{1};
        
        trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
        [~, trial_frames] = f_dv_compute_window_t(trial_window, mean(cat(1,cdata.volume_period)));
        
        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_frames);

        [selected_cells, ~, ~, ~, resp_cells] = f_dv_get_resp_vals_cells(app, stats1, tn_all);

        for n_tn = 1:num_tn
            tn1 = tn_all(n_tn);
            
            tr_idx = logical(sum(trial_types == app.ops.context_types_all(tn1)',2));
            tr_data = trial_data_sort(:,:,tr_idx);
            
            num_trials = sum(tr_idx);

            selected_cells2 = selected_cells(:,n_tn);
            resp_cells2 = resp_cells(:,n_tn);
            
            if sum(resp_cells2)>5
                
                tr_data2 = tr_data(selected_cells2,:,:);
                %firing_rate2 = firing_rate(resp_cells,:);
        
                hc_params.plot_dist_mat = 0;
                hc_params.plot_clusters = 0;
                hc_params.num_clust = 1;
        
        
                tr_data_2d_tr = reshape(tr_data2, [], num_trials);
                hclust_out_trial = f_hcluster_wrap(tr_data_2d_tr', hc_params);
                %tr_data3 = tr_data2(:,:,hclust_out_trial.dend_order);
        
                SI = 1-hclust_out_trial.dist;
                SI_vals = tril(SI,-1);
                SI_vals(SI_vals==0) = [];
                
                corr_vals(n_dset, n_tn) = mean(SI_vals);
            end
       
        end
    end
end


isi_uq = unique(isi_vals);
num_isi = numel(isi_uq);

color1 = jet(10);
figure; hold on
for n_freq = 1:num_tn
    corr_temp = zeros(num_isi,1);
    for n_isi = 1:num_isi
        isi_idx = isi_vals == isi_uq(n_isi);
        temp_data = corr_vals(isi_idx,n_freq);
        corr_temp(n_isi) = mean(temp_data(~isnan(temp_data)));
    end
    plot(isi_uq, corr_temp, 'o-', 'color', color1(tn_all(n_freq),:), 'linewidth', 2)
end
xlabel('ISI duration'); ylabel('Pairwise correlation')
title(sprintf('Mean pairwise correlations; %s', title_tag), 'interpreter', 'none');
xlim([0, 4.5]);


corr_mean = zeros(num_isi,1);
corr_sem = zeros(num_isi,1);
indiv_dat = cell(num_isi,1);
lab_all = cell(num_isi,1);
lab_all2 = cell(num_isi,1);
for n_isi = 1:num_isi
    isi_idx = isi_vals == isi_uq(n_isi);
    temp_data = corr_vals(isi_idx,:);
    temp_data2 = temp_data(~isnan(temp_data));
    indiv_dat{n_isi} = temp_data2;
    corr_mean(n_isi) = mean(temp_data2);
    corr_sem(n_isi) = std(temp_data2)/sqrt(numel(temp_data2)-1);
    lab_all{n_isi} = repmat({num2str(round(isi_uq(n_isi),1))}, [numel(temp_data2), 1]);
    lab_all2{n_isi} = num2str(round(isi_uq(n_isi),1));
end

figure; hold on
for n_isi = 1:num_isi
    x_data = (rand(numel(indiv_dat{n_isi}),1)-0.5)/8 + isi_uq(n_isi);
    %x_data = zeros(numel(indiv_dat{n_isi}),1) + n_isi;
    plot(x_data, indiv_dat{n_isi}, '.', color=[0.4 0.4 0.4]);
end
%plot(isi_uq, corr_mean, 'o-', 'linewidth', 2)
errorbar(isi_uq, corr_mean, corr_sem, '.-', 'linewidth', 2, color='k', markersize=15)
xlabel('ISI duration'); 
ylabel('Pairwise correlation')
title(sprintf('Mean pairwise correlations; %s; %s', title_tag, app.ResponsivecellsselectDropDown.Value), 'interpreter', 'none');
xlim([0, 4.5]);
ylim([0 0.5]);

figure; hold on
pl2 = cell(2,1);
for n_isi = 1:num_isi
    x_data = (rand(numel(indiv_dat{n_isi}),1)-0.5)/8 + isi_uq(n_isi);
    %x_data = zeros(numel(indiv_dat{n_isi}),1) + n_isi;
    plot(x_data, indiv_dat{n_isi}, '.', color=[0.4 0.4 0.4]);
end
pl2{1} = errorbar(isi_uq, corr_mean, corr_sem, '.', 'linewidth', 2, color='k', markersize=15);
fit_data = f_get_fit(isi_uq, corr_mean, 'exp');

x_fit = linspace(isi_uq(1), isi_uq(end), 100);
y_fit = fit_data.fit_eq(x_fit, fit_data.fit_pars);
pl2{2} = plot(x_fit, y_fit, '--k');
legend([pl2{:}], {'data', 'exp fit'});
xlabel('ISI duration'); 
ylabel('Pairwise correlation')
title(sprintf('Mean pairwise correlations; %s; %s; tau=%.3f, %.3fstd', title_tag, app.ResponsivecellsselectDropDown.Value, fit_data.fit_pars(2), sqrt(fit_data.param_cov(2,2))), 'interpreter', 'none');
xlim([0, 4.5]);
ylim([0 0.5]);
[p_all, tbl_all, stats_all]  = anova1(cat(1, indiv_dat{:}), cat(1,lab_all{:}), 'off');

title_tag5 = sprintf('between isi pts; %s', title_tag); 

f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag5, lab_all2, 1, 0);

% [h,p,ci,stats] = ttest2(indiv_dat{1}, cat(1,indiv_dat{2:4}))


%figure; imagesc(reshape(color1, [10 1 3]))

end