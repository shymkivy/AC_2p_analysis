function f_dv_decoder_onevone(app)

decoder_type = 'bayes'; % tree, svm, bayes
by_frame = 1;

tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';
[num_gr, num_tn] = size(tn_all);

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = size(data,1);
params = f_dv_gather_params(app);

[region_num, reg_tag] = f_dv_get_region_sel_val(app);
num_regions = numel(region_num);
reg_all = app.ops.regions_to_analyze;

ddata = data(1,:);
[cdata, ~] = f_dv_get_new_cdata_stats(app, ddata, params);

if by_frame
    trial_window = [-1, 3];
    [plot_t, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1).volume_period);
    dec_acc_frames = zeros(num_dsets,sum(trial_frames), num_gr, num_regions);
    dec_acc_frames_shuff = zeros(num_dsets,sum(trial_frames), num_gr, num_regions);
    dec_acc_frames_bycl = zeros(num_dsets,sum(trial_frames), num_gr, num_regions, num_tn);
    dec_acc_frames_bycl_shuff = zeros(num_dsets,sum(trial_frames), num_gr, num_regions, num_tn);
else
    trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
    [plot_t, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1).volume_period);
    dec_acc_curr_last_shuff = zeros(num_dsets,3, num_gr, num_regions);
end

title_tag2 = sprintf('%s; %s reg; %dms smooth', title_tag, reg_tag, app.SmoothsigmamsEditField.Value);

mouse_id = cell(num_gr, num_dsets);
dset_id = zeros(num_gr, num_dsets);
num_cells_avail = zeros(num_gr, num_dsets, num_regions);
done_dec = false(num_gr, num_dsets, num_regions);

fprintf('dset n/%d: ', num_dsets);

for n_dset = 1:num_dsets
    
    fprintf('%d..', n_dset);

    ddata = data(n_dset,:);
    [cdata, stats1] = f_dv_get_new_cdata_stats(app, ddata, params);
    mmn_freq = ddata.MMN_freq{1};
    stim_times = ddata.stim_frame_index{1};

    firing_rate = cat(1,cdata.S_sm);
    
    trial_types = ddata.trial_types{1};
    

    num_cells = sum([stats1.num_cells]);

    if app.UseregdatalabelsCheckBox.Value
        if ~isempty(ddata.registered_data{1})
            reg_cell_labels = ddata.registered_data{1}.reg_labels;
        else
            reg_cell_labels = zeros(num_cells,1);
        end
    else
        reg_idx = find(strcmpi(reg_all, ddata.area));
        reg_cell_labels = ones(num_cells,1)*reg_idx;
    end

    for n_gr = 1:num_gr
        tn1 = tn_all(n_gr, :);
        mouse_id{n_gr, n_dset} = ddata.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;

        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn1);
        resp_cells2 = logical(sum(resp_cells,2));
        for n_reg = 1:num_regions
            reg_idx = reg_cell_labels == region_num(n_reg);
            resp_reg_cell = and(resp_cells2, reg_idx);
            
            num_cells_avail(n_gr, n_dset, n_reg) = sum(resp_reg_cell);
            
            if num_cells_avail(n_gr, n_dset, n_reg) > 10

                done_dec(n_gr, n_dset, n_reg) = 1;
    
                firing_rate2 = firing_rate(resp_reg_cell,:);
    
                trial_data_sort = f_get_stim_trig_resp(firing_rate2, stim_times, trial_frames);
    
                trial_types2 = trial_types(2:400);
                trial_types2_last = trial_types(1:399);
                trial_types2_shuff = trial_types2(randperm(numel(trial_types2)));
                
                trial_data_sort2 = trial_data_sort(:,:,2:400);
    
                if by_frame
                    for n_fr = 1:sum(sum(trial_frames))
                        
                        trial_data_sort3 = squeeze(trial_data_sort2(:,n_fr,:));
                        
                        trial_data_sort3n = trial_data_sort3 - mean(trial_data_sort3(:));
                        trial_data_sort3n = trial_data_sort3n/std(trial_data_sort3n(:));
                    
                        if strcmpi(decoder_type, 'tree')
                            [trainedClassifier_tree, dec_acc_frames(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl(n_dset,n_fr,n_gr,n_reg,:)] = decoder_tree(trial_data_sort3n', trial_types2);
                            [trainedClassifier_tree_shuff, dec_acc_frames_shuff(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl_shuff(n_dset,n_fr,n_gr,n_reg,:)] = decoder_tree(trial_data_sort3n', trial_types2_shuff);
                        elseif strcmpi(decoder_type, 'bayes')
                            [trainedClassifier_bayes, dec_acc_frames(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl(n_dset,n_fr,n_gr,n_reg,:)] = decoder_naivebayes(trial_data_sort3n', trial_types2);
                            [trainedClassifier_bayes_shuff, dec_acc_frames_shuff(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl_shuff(n_dset,n_fr,n_gr,n_reg,:)] = decoder_naivebayes(trial_data_sort3n', trial_types2_shuff);
                        elseif strcmpi(decoder_type, 'svm')
                            [trainedClassifier_svm, dec_acc_frames(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl(n_dset,n_fr,n_gr,n_reg,:)] = decoder_svm(trial_data_sort3n', trial_types2, 1);
                            [trainedClassifier_svm_shuff, dec_acc_frames_shuff(n_dset,n_fr,n_gr,n_reg), dec_acc_frames_bycl_shuff(n_dset,n_fr,n_gr,n_reg,:)] = decoder_svm(trial_data_sort3n', trial_types2_shuff, 1);
                        end
                    end
                else
                    trial_data_sort3 = squeeze(mean(trial_data_sort2,2));
               
                    trial_data_sort3n = trial_data_sort3 - mean(trial_data_sort3(:));
                    trial_data_sort3n = trial_data_sort3n/std(trial_data_sort3n(:));
                
                    if sprintf(decoder_type, 'tree')
                        [trainedClassifier_tree, dec_acc_curr_last_shuff(n_dset,1,n_gr,n_reg)] = decoder_tree(trial_data_sort3n', trial_types2);
                        [trainedClassifier_tree_last, dec_acc_curr_last_shuff(n_dset,2,n_gr,n_reg)] = decoder_tree(trial_data_sort3n', trial_types2_last);
                        [trainedClassifier_tree_shuff, dec_acc_curr_last_shuff(n_dset,3,n_gr,n_reg)] = decoder_tree(trial_data_sort3n', trial_types2_shuff);
                    elseif sprintf(decoder_type, 'bayes')
                        [trainedClassifier_bayes, dec_acc_curr_last_shuff(n_dset,1,n_gr,n_reg)] = decoder_naivebayes(trial_data_sort3n', trial_types2);
                        [trainedClassifier_bayes_last, dec_acc_curr_last_shuff(n_dset,2,n_gr,n_reg)] = decoder_naivebayes(trial_data_sort3n', trial_types2_last);
                        [trainedClassifier_bayes_shuff, dec_acc_curr_last_shuff(n_dset,3,n_gr,n_reg)] = decoder_naivebayes(trial_data_sort3n', trial_types2_shuff);
                    elseif sprintf(decoder_type, 'svm')
                        [trainedClassifier_svm, dec_acc_curr_last_shuff(n_dset,1,n_gr,n_reg)] = decoder_svm(trial_data_sort3n', trial_types2, 1);
                        [trainedClassifier_svm_last, dec_acc_curr_last_shuff(n_dset,2,n_gr,n_reg)] = decoder_svm(trial_data_sort3n', trial_types2_last, 1);
                        [trainedClassifier_svm_shuff, dec_acc_curr_last_shuff(n_dset,3,n_gr,n_reg)] = decoder_svm(trial_data_sort3n', trial_types2_shuff, 1);
                    end
                end
                % % mv regression
                % MnrModel_cont = fitmnr(trial_data_sort3n', trial_types2);
                % 
                % 
                % [d,p,stats] = manova1(trial_data_sort3n', trial_types2);
                % 
                % load fisheriris
                % 
                % MnrModel = fitmnr(meas,species);
                % MnrModel.Coefficients
            
                % beta = mvregress(X,Y)
            end
        end
    end  
end

fprintf('\n');



for n_gr = 1:num_gr
    for n_reg = 1:num_regions
        done_dec2 = done_dec(n_gr, :, n_reg);
        dec_acc_frames2 = dec_acc_frames(done_dec2,:,n_gr,n_reg);
        dec_acc_frames_shuff2 = dec_acc_frames_shuff(done_dec2,:,n_gr,n_reg);
        if by_frame
            if ~isempty(dec_acc_frames_shuff2)
                figure; hold on; axis tight
                plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
                plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2)
            
                plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
                plot(plot_t, mean(dec_acc_frames2,1), color=[0 0.4470 0.7410], LineWidth=2)
                title(sprintf('%s; %s decoder, freqs; %s', reg_all{region_num(n_reg)}, decoder_type, title_tag2), 'interpreter', 'none')
            end
        else
            mean(dec_acc_curr_last_shuff)
        end
    end
end

for n_gr = 1:num_gr
    figure; hold on; axis tight
    pl_all = cell(num_regions+1,1);
    has_reg_data = false(num_regions+1,1);
    for n_reg = 1:num_regions
        done_dec2 = done_dec(n_gr, :, n_reg);
        dec_acc_frames2 = dec_acc_frames(done_dec2,:,n_gr,n_reg);
        dec_acc_frames_shuff2 = dec_acc_frames_shuff(done_dec2,:,n_gr,n_reg);
        if by_frame
            if ~isempty(dec_acc_frames_shuff2)
                has_reg_data(n_reg) = 1;
                has_reg_data(num_regions+1) = 1;
                %plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
                pl_all{num_regions+1} = plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2);
                col2 = app.ops.cond_colors{region_num(n_reg)};
                %plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
                pl_all{n_reg} = plot(plot_t, mean(dec_acc_frames2,1), color=col2, LineWidth=2);
            end
        else
            mean(dec_acc_curr_last_shuff)
        end
    end
    title(sprintf('%s decoder, freqs; %s', decoder_type, title_tag2), 'interpreter', 'none')
    legend([pl_all{has_reg_data}], [reg_all(region_num(has_reg_data(1:num_regions))); {'Shuffle'}])
end

for n_gr = 1:num_gr
    for n_reg = 1:num_regions
        done_dec2 = done_dec(n_gr, :, n_reg);
        dec_acc_frames_bycl2 = dec_acc_frames_bycl(done_dec2,:,n_gr,n_reg,:);
        dec_acc_frames_bycl_shuff2 = dec_acc_frames_bycl_shuff(done_dec2,:,n_gr,n_reg,:);
        if by_frame
            if ~isempty(dec_acc_frames_bycl2)
                figure; hold on; axis tight
                pl_all = cell(num_tn+1);
                for n_tt = 1:num_tn
                    %plot(plot_t, dec_acc_frames_bycl_shuff2', color=[0 0 0 0.2])
                    pl_all{num_tn+1} = plot(plot_t, mean(dec_acc_frames_bycl_shuff2(:,:,:,:,n_tt),1), color=[0 0 0], LineWidth=2);
                end
                for n_tt = 1:num_tn
                    col2 = app.ops.context_types_all_colors2{tn_all(n_tt)};
                    %plot(plot_t, dec_acc_frames_bycl2', color=[0 0.4470 0.7410 0.2])
                    pl_all{n_tt} = plot(plot_t, mean(dec_acc_frames_bycl2(:,:,:,:,n_tt),1), color=col2, LineWidth=2);
                end
                title(sprintf('%s; %s decoder, freqs; %s', reg_all{region_num(n_reg)}, decoder_type, title_tag2), 'interpreter', 'none');
                legend([pl_all{:}], [app.ops.context_types_labels(tn_all); {'Shuffle'}]);
            end
        else
            mean(dec_acc_curr_last_shuff)
        end
    end
end

end