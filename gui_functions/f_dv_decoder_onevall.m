function f_dv_decoder_onevall(app)

decoder_type = 'svm'; % tree, svm, bayes

trial_num_selection = 'min'; % all, median, mean, min

tn_all = f_dv_get_trial_number(app);
%tt_all = app.ops.context_types_all(tn_all)';
[num_gr, num_tn] = size(tn_all);

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = size(data,1);
params = f_dv_gather_params(app);

[region_num, reg_tag] = f_dv_get_region_sel_val(app);
num_regions = numel(region_num);
reg_all = app.ops.regions_to_analyze;

ddata = data(1,:);
[cdata, ~] = f_dv_get_new_cdata_stats(app, ddata, params);

trial_window = [-1, 3];
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, cdata(1).volume_period);
dec_acc_frames = zeros(num_dsets,sum(trial_frames), num_gr, num_regions);
dec_acc_frames_shuff = zeros(num_dsets,sum(trial_frames), num_gr, num_regions);
dec_acc_frames_bycl = zeros(num_dsets,sum(trial_frames), num_gr, num_regions, num_tn);
dec_acc_frames_bycl_shuff = zeros(num_dsets,sum(trial_frames), num_gr, num_regions, num_tn);


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
    if ~isempty(mmn_freq)
        trial_types_ctx2 = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);
        trial_types_all = [trial_types, trial_types_ctx2];
    else
        trial_types_all = trial_types;
    end

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
                tt1 = reshape(app.ops.context_types_all(tn1), [1, 1, numel(tn1)]);

                idx1 = logical(sum(trial_types_all == tt1,3));

                num_tt2 = size(trial_types_all,2);
                trial_types2c = cell(num_tt2,1);
                trial_data_sort2c = cell(num_tt2,1);
                for n_r = 1:num_tt2
                    trial_types2c{n_r} = trial_types_all(idx1(:,n_r),n_r);
                    trial_data_sort2c{n_r} = trial_data_sort(:,:,idx1(:,n_r));
                end

                trial_types2 = cat(1, trial_types2c{:});
                trial_data_sort2 = cat(3,trial_data_sort2c{:});
                
                trial_counts = sum(trial_types2 == app.ops.context_types_all(tn1)',1);
                if strcmpi(trial_num_selection, 'all')
                    max_tr = max(trial_counts);
                elseif strcmpi(trial_num_selection, 'mean')
                    max_tr = round(mean(trial_counts));
                elseif strcmpi(trial_num_selection, 'median')
                    max_tr = round(defian(trial_counts));
                elseif strcmpi(trial_num_selection, 'min')
                    max_tr = round(min(trial_counts));
                end

                trial_types3c = cell(num_tn,1);
                trial_data_sort3c = cell(num_tn,1);
                for n_tt = 1:num_tn
                    tt2 = app.ops.context_types_all(tn1(n_tt));
                    tt_idx1 = trial_types2 == tt2;
                    if sum(tt_idx1) > max_tr
                        samp_tr1 = randsample(find(tt_idx1), max_tr);
                        trial_types3c{n_tt} = trial_types2(samp_tr1);
                        trial_data_sort3c{n_tt} = trial_data_sort2(:,:,samp_tr1);
                    else
                        trial_types3c{n_tt} = trial_types2(tt_idx1);
                        trial_data_sort3c{n_tt} = trial_data_sort2(:,:,tt_idx1);
                    end
                end

                trial_types3 = cat(1, trial_types3c{:});
                trial_data_sort3 = cat(3,trial_data_sort3c{:});
                
                sh_idx1 = randperm(numel(trial_types3));
                
                trial_types4 = trial_types3(sh_idx1);
                trial_data_sort4 = trial_data_sort3(:,:,sh_idx1);

                trial_types4_shuff = trial_types4(randperm(numel(trial_types4)));
                
                for n_fr = 1:sum(sum(trial_frames))
                    
                    trial_data_sort4n = squeeze(trial_data_sort4(:,n_fr,:));
                    trial_data_sort4n = trial_data_sort4n - mean(trial_data_sort4n(:));
                    trial_data_sort4n = trial_data_sort4n/std(trial_data_sort4n(:));
                
                    if strcmpi(decoder_type, 'tree')
                        dec_out = decoder_tree(trial_data_sort4n', trial_types4);
                        dec_out_shuff = decoder_tree(trial_data_sort4n', trial_types4_shuff);
                    elseif strcmpi(decoder_type, 'bayes')
                        dec_out = decoder_naivebayes(trial_data_sort4n', trial_types4);
                        dec_out_shuff = decoder_naivebayes(trial_data_sort4n', trial_types4_shuff);
                    elseif strcmpi(decoder_type, 'svm')
                        dec_out = decoder_svm(trial_data_sort4n', trial_types4, 1);
                        dec_out_shuff = decoder_svm(trial_data_sort4n', trial_types4_shuff, 1);
                    end

                    dec_acc_frames(n_dset,n_fr,n_gr,n_reg) = dec_out.validationAccuracy;
                    dec_acc_frames_bycl(n_dset,n_fr,n_gr,n_reg,:) = dec_out.acc_by_class;
                    dec_acc_frames_shuff(n_dset,n_fr,n_gr,n_reg) = dec_out_shuff.validationAccuracy;
                    dec_acc_frames_bycl_shuff(n_dset,n_fr,n_gr,n_reg,:) = dec_out_shuff.acc_by_class;
                end
                
            end
        end
    end  
end

fprintf('\n done \n');


for n_gr = 1:num_gr
    for n_reg = 1:num_regions
        done_dec2 = done_dec(n_gr, :, n_reg);
        dec_acc_frames2 = dec_acc_frames(done_dec2,:,n_gr,n_reg);
        dec_acc_frames_shuff2 = dec_acc_frames_shuff(done_dec2,:,n_gr,n_reg);

        if ~isempty(dec_acc_frames_shuff2)
            figure; hold on; axis tight
            plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
            plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2)
        
            plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
            plot(plot_t, mean(dec_acc_frames2,1), color=[0 0.4470 0.7410], LineWidth=2)
            title(sprintf('%s; %s decoder, freqs; %s', reg_all{region_num(n_reg)}, decoder_type, title_tag2), 'interpreter', 'none')
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

        if ~isempty(dec_acc_frames_shuff2)
            has_reg_data(n_reg) = 1;
            has_reg_data(num_regions+1) = 1;
            %plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
            pl_all{num_regions+1} = plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2);
            col2 = app.ops.cond_colors{region_num(n_reg)};
            %plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
            pl_all{n_reg} = plot(plot_t, mean(dec_acc_frames2,1), color=col2, LineWidth=2);
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

        if ~isempty(dec_acc_frames_bycl2)
            figure; hold on; axis tight
            pl_all = cell(num_tn+1);
            for n_tt = 1:num_tn
                %plot(plot_t, dec_acc_frames_bycl_shuff2', color=[0 0 0 0.2])
                pl_all{num_tn+1} = plot(plot_t, mean(dec_acc_frames_bycl_shuff2(:,:,:,:,n_tt),1), color=[.5 .5 .5], LineWidth=2);
            end
            for n_tt = 1:num_tn
                col2 = app.ops.context_types_all_colors2{tn_all(n_tt)};
                %plot(plot_t, dec_acc_frames_bycl2', color=[0 0.4470 0.7410 0.2])
                pl_all{n_tt} = plot(plot_t, mean(dec_acc_frames_bycl2(:,:,:,:,n_tt),1), color=col2, LineWidth=2);
            end
            title(sprintf('%s; %s decoder, freqs; %s', reg_all{region_num(n_reg)}, decoder_type, title_tag2), 'interpreter', 'none');
            legend([pl_all{:}], [app.ops.context_types_labels(tn_all); {'Shuffle'}]);
        end
    end
end

end