function f_dv_plot_single_trial_PCA(app)

normalize1 = 0;
plot_extra = app.plotsuperdeetsCheckBox.Value;

num_comp = app.plotPCAdimSpinner.Value;
num_pl_d = str2double(app.plotaxesDropDown.Value);
num_pl = floor(num_comp/num_pl_d);
num_plot_comp = num_pl*num_pl_d;

trs_all = reshape(1:num_plot_comp, num_pl_d, []);
%%
n_pl = app.mplSpinner.Value;

ddata = app.ddata;
title_tag = ['plane ' app.ddata.dset_name_full{1}];
num_dsets = 1;

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);
num_t = sum(trial_frames);

tn0 = f_dv_get_trial_number(app, ddata.MMN_freq{1});

[num_tn_gr, num_tn] = size(tn0);

params = f_dv_gather_params(app);

%reg_all = app.ops.regions_to_analyze;
[region_num_all, reg_tag, leg_list] = f_dv_get_region_sel_val2(app);

num_reg = size(region_num_all,1);

mouse_id = cell(num_tn_gr, num_dsets);
dset_id = zeros(num_tn_gr, num_dsets);

resp_all_mean = cell(num_tn_gr, num_dsets, num_reg);
resp_all_trials = cell(num_tn_gr, num_dsets, num_reg, num_tn);
resp_num_trials = zeros(num_tn_gr, num_dsets, num_reg, num_tn);
cell_counts = zeros(num_tn_gr, num_dsets, num_reg);

fprintf('dset #/%d: ', num_dsets);
for n_dset = 1:num_dsets
    fprintf('..%d', n_dset);
    data1 =  ddata(n_dset,:);
    stats1 = data1.stats{n_pl};
    params.n_dset = find(data1.idx == app.data.idx);
    tn1 = f_dv_get_trial_number(app, data1.MMN_freq{1});

    cdata = f_dv_compute_cdata(data1, params);

    firing_rate = cdata.S_sm;
    mmn_freq = data1.MMN_freq{1};
    trial_types = data1.trial_types{1};
    trial_types_ctx = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);

    stim_times = data1.stim_frame_index{n_pl};
    
    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);

    reg_cell_labels = f_dv_get_area_label(app, data1);
    
    for n_tngr = 1:num_tn_gr
        mouse_id{n_tngr, n_dset} = data1.mouse_id{1};
        dset_id(n_tngr, n_dset) = n_dset;

        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn1(n_tngr,:));
        %cell_is_resp = stats1.peak_resp_cells(:,tn_all);

        for n_reg = 1:num_reg
            region_num = region_num_all(n_reg,:);

            reg_cell_idx = logical(sum(reg_cell_labels == region_num,2));

            % get resp cells
            resp_cell_idx = logical(sum(resp_cells,2).*reg_cell_idx);
            num_cells = sum(resp_cell_idx);

            if num_cells
                cell_counts(n_tngr, n_dset, n_reg) = num_cells;
                resp_all_mean{n_tngr, n_dset, n_reg} = zeros(num_cells, num_t, num_tn);
                for n_tn = 1:num_tn
                    idx1 = logical(sum([trial_types, trial_types_ctx] == app.ops.context_types_all(tn1(n_tn)),2));
                    temp_resp = trial_data_sort(resp_cell_idx,:,idx1);
                    resp_all_mean{n_tngr, n_dset, n_reg}(:,:,n_tn) = mean(temp_resp,3);
                    resp_all_trials{n_tngr, n_dset, n_reg, n_tn} = temp_resp;
                    resp_num_trials(n_tngr, n_dset, n_reg, n_tn) = size(temp_resp,3);
                end
            end
        end
    end
end
fprintf('\n');

num_effdsets = num_dsets*num_tn_gr;

resp_all2 = reshape(resp_all_mean, num_effdsets, num_reg);
resp_all_trials2 = reshape(resp_all_trials, num_effdsets, num_reg, num_tn);
resp_num_trials2 = reshape(resp_num_trials, num_effdsets, num_reg, num_tn);
cell_counts2  = reshape(cell_counts, num_effdsets, num_reg);
mouse_id2 = reshape(mouse_id, num_effdsets, 1);
dset_id2 = reshape(dset_id, num_effdsets, 1);

[gr_id, gr_tag] = f_dv_combine_data(app, mouse_id2, dset_id2);

%% now the pca

gr_all = unique(gr_id);
num_gr = numel(gr_all);

for n_dset = 1:num_dsets
    top_comp_comb = cell(1, num_reg);
    exp_var_comb = cell(1, num_reg);
    for n_reg = 1:num_reg
        num_cells = sum(cell_counts2(:,n_reg));
        if num_cells
            
            resp2_tr_mean = cat(3,resp_all2{n_dset,n_reg,:});
            
            [num_cells1, num_t1, num_tr1] = size(resp2_tr_mean);

            resp2d_tr_mean = reshape(resp2_tr_mean, num_cells, num_t1*num_tr1);
            
            [~,score,~,~,explained,mu] = pca(resp2d_tr_mean');

            resp2 = cat(3,resp_all_trials2{n_dset,n_reg,:});

            resp_num_trials3 = squeeze(resp_num_trials2(n_dset,n_reg,:));
            
            resp2d = reshape(resp2, num_cells, []);
            
            if normalize1
                resp2d = f_normalize(resp2d, 'norm_mean_std');
            end
        
            %[U,S,V] = svd(resp_all2dn);
            % score*coeff'
            %warning('off', 'stats:pca:ColRankDefX')
            [~,score,~,~,explained,mu] = pca(resp2d');
            
            top_comp = score(:,1:num_comp);
            top_comp2 = reshape(top_comp, [num_t, sum(resp_num_trials3), num_comp]);
            
            top_comp_comb{1, n_reg} = top_comp2;
            exp_var_comb{1, n_reg} = explained;
        end
    end
    
    colors_tn = app.ops.context_types_all_colors2;
    %tn1 = tn_all(1, :);
    
    for n_reg = 1:num_reg
        num_cells = sum(cell_counts2(:,n_reg));
        if num_cells
            resp_num_trials3 = squeeze(resp_num_trials2(n_dset,n_reg,:));
    
            top_comp2 = top_comp_comb{1, n_reg};
            exp_var2 = exp_var_comb{1, n_reg};
            
            num_tr = sum(resp_num_trials3);
            tr_idx1 = cell(num_tn, 1);
            for n_tn = 1:num_tn
                if iscell(tn1)
                    tn1 = tn1{n_dset};
                else
                    tn1 = tn1;
                end
                tr_idx1{n_tn} = ones(resp_num_trials3(n_tn), 1) * tn1(n_tn);
            end
            tr_idx2 = cat(1, tr_idx1{:});
            
            figure;
            plot(exp_var2, 'ko-')

            for n_pl = 1:num_pl
                trs1 = trs_all(:, n_pl);
                if num_pl_d == 2
                    figure; hold on;
                    for n_tr = 1:num_tr
                        plot(top_comp2(:,n_tr,trs1(1)), top_comp2(:,n_tr,trs1(2)), 'color', colors_tn{tr_idx2(n_tr)})
                    end
                    xlabel(sprintf('PC %d', trs1(1)))
                    ylabel(sprintf('PC %d', trs1(2)))
                
                elseif num_pl_d == 3
                    figure; hold on;
                    for n_tr = 1:num_tr
                        plot3(top_comp2(:,n_tr,trs1(1)), top_comp2(:,n_tr,trs1(2)), top_comp2(:,n_tr,trs1(3)), 'color', colors_tn{tr_idx2(n_tr)})
                    end
                    xlabel(sprintf('PC %d', trs1(1)))
                    ylabel(sprintf('PC %d', trs1(2)))
                    zlabel(sprintf('PC %d', trs1(3)))
                end
            end

        end
    end
end
end