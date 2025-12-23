function f_plot_decoder_data_by_regions(dec_data, plot_t, region_id, reg_labels, trial_group_id, params, ops)

dec_acc_frames = cat(1,dec_data.accuracy);
dec_acc_frames_shuff = cat(1,dec_data.accuracy_shuff);

dec_acc_frames_bycl = cat(1,dec_data.accuracy_by_class);
dec_acc_frames_bycl_shuff = cat(1,dec_data.accuracy_by_class_shuff);

regions = unique(region_id);
num_regions = numel(regions);
groups = unique(trial_group_id);
num_gr = numel(groups);

tn_all = f_dv_get_trial_number(params);

for n_gr = 1:num_gr
    for n_reg = 1:num_regions
        dec_idx = and(region_id==n_reg, trial_group_id == n_gr);
        dec_acc_frames2 = dec_acc_frames(dec_idx,:);
        dec_acc_frames_shuff2 = dec_acc_frames_shuff(dec_idx,:);
        if ~isempty(dec_acc_frames_shuff2)
            figure; hold on; axis tight
            plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
            plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2)
        
            plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
            plot(plot_t, mean(dec_acc_frames2,1), color=[0 0.4470 0.7410], LineWidth=2)
            title(sprintf('%s; %s decoder, freqs', reg_labels{n_reg}, params.decoder_type), 'interpreter', 'none');
            xlim([-0.5, 2.5]);
        end
    end
end

if num_regions > 1
    colors2 = ops.cond_colors;
else
    colors2 = {[0 0.4470 0.7410]};
end

for n_gr = 1:num_gr
    figure; hold on; axis tight
    pl_all = cell(num_regions+1,1);
    has_reg_data = false(num_regions+1,1);
    for n_reg = 1:num_regions
        dec_idx = and(region_id==n_reg, trial_group_id == n_gr);
        dec_acc_frames2 = dec_acc_frames(dec_idx,:);
        dec_acc_frames_shuff2 = dec_acc_frames_shuff(dec_idx,:);

        if ~isempty(dec_acc_frames_shuff2)
            has_reg_data(n_reg) = 1;
            has_reg_data(num_regions+1) = 1;
            %plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
            pl_all{num_regions+1} = plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2);
            col2 = colors2{n_reg};
            %plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
            pl_all{n_reg} = plot(plot_t, mean(dec_acc_frames2,1), color=col2, LineWidth=2);
        end
    end
    title(sprintf('%s decoder, freqs', params.decoder_type), 'interpreter', 'none');
    legend([pl_all{has_reg_data}], [reg_labels(has_reg_data(1:num_regions)), {'Shuffle'}]);
    xlim([-0.5, 2.5]);
end

sig_plot = [0.001, 0.01, 0.05];
sig_range = [0.75 0.95];

sig_space = diff(sig_range)/(numel(sig_plot)+1);
reg_space = sig_space/num_regions/2;

for n_gr = 1:num_gr
    figure; hold on; axis tight
    pl_all = cell(num_regions+1,1);
    has_reg_data = false(num_regions+1,1);

    if params.plot_stim
        for n_st = 1:3
            r1 = rectangle('Position', [n_st-1 0 0.5 1]);
            if isprop(r1, "FaceAlpha")
                r1.FaceColor = [ops.context_types_all_colors2{params.stim_freq_color}];
                r1.FaceAlpha = params.stim_transparancy;
            else
                r1.FaceColor = [ops.context_types_all_colors2{params.stim_freq_color} params.stim_transparancy];
            end
            r1.EdgeColor = [ops.context_types_all_colors2{params.stim_freq_color} params.stim_transparancy];
        end
    end

    shuff_all = cell(num_regions,1);
    for n_reg = 1:num_regions
        %done_dec2 = done_dec(n_gr, :, n_reg);
        dec_idx = and(region_id==n_reg, trial_group_id == n_gr);
        dec_acc_frames2 = dec_acc_frames(dec_idx,:);
        dec_acc_frames_shuff2 = dec_acc_frames_shuff(dec_idx,:);
        shuff_all{n_reg} = dec_acc_frames_shuff2;

        if ~isempty(dec_acc_frames_shuff2)
            has_reg_data(n_reg) = 1;
            has_reg_data(num_regions+1) = 1;
            %plot(plot_t, dec_acc_frames_shuff2', color=[0 0 0 0.2])
            %pl_all{num_regions+1} = plot(plot_t, mean(dec_acc_frames_shuff2,1), color=[0 0 0], LineWidth=2);
            col2 = colors2{n_reg};
            %plot(plot_t, dec_acc_frames2', color=[0 0.4470 0.7410 0.2])
            num_dec1 = size(dec_acc_frames2,1);
            if num_dec1 > 1
                s1 = shadedErrorBar_YS(plot_t, mean(dec_acc_frames2,1), std(dec_acc_frames2,[],1)./sqrt(num_dec1-1), col2);
                pl_all{n_reg} = s1.mainLine;
            else
                pl_all{n_reg} = plot(plot_t, mean(dec_acc_frames2,1), color=col2);
            end
            %pl_all{n_reg} = plot(plot_t, mean(dec_acc_frames2,1), color=col2, LineWidth=2);
        end
    end
    shuff_all2 = cat(1, shuff_all{:});
    num_dec2 = size(shuff_all2,1);
    s1 = shadedErrorBar_YS(plot_t, mean(shuff_all2,1), std(shuff_all2,[],1)./sqrt(num_dec2-1), [0 0 0]);
    pl_all{num_regions+1} = s1.mainLine;

    for n_reg = 1:num_regions
        %done_dec2 = done_dec(n_gr, :, n_reg);
        dec_idx = and(region_id==n_reg, trial_group_id == n_gr);
        dec_acc_frames2 = dec_acc_frames(dec_idx,:);
        
        samp11 = dec_acc_frames2;
        samp22 = shuff_all2;

        %samp1 = dec_acc_frames2(:,1);
        %samp2 = shuff_all2(:,1);
        %[h,p,ci,stats] = ttest2(samp1, samp2)
        
        n1 = size(samp11,1);
        n2 = size(samp22,1);
        
        t_vals1 = (mean(samp11, 1) - mean(samp22, 1))./sqrt(var(samp11, [] ,1)/n1 + var(samp22, [], 1)/n2);
        df1 = n1 + n2 - 2;
        p_vals1 = (1 - tcdf(abs(t_vals1), df1))*2;
        
        col2 = colors2{n_reg};
 
        for n_sig = 1:numel(sig_plot)
            idx1 = p_vals1 < sig_plot(n_sig);
            sig_trace = nan(numel(plot_t),1);
            sig_trace(idx1) = 1;

            y_level = max(sig_range) - (n_sig-1)*sig_space - (n_reg-1)*reg_space;
            plot(plot_t, sig_trace*y_level, '.-', color=col2)

            if n_reg == 1
                text(-0.5+sig_space, y_level+reg_space, ['p<' num2str(sig_plot(n_sig))]);
            end
        end

    end
    legend([pl_all{has_reg_data}], [reg_labels(has_reg_data(1:num_regions)), {'Shuffle'}]);
    title(sprintf('%s decoder, freqs', params.decoder_type), 'interpreter', 'none');
    xlim([-0.5, 2.5]);
end

num_tn = size(dec_acc_frames_bycl,3);
for n_gr = 1:num_gr
    for n_reg = 1:num_regions
        %done_dec2 = done_dec(n_gr, :, n_reg);
        dec_idx = and(region_id==n_reg, trial_group_id == n_gr);
        dec_acc_frames_bycl2 = dec_acc_frames_bycl(dec_idx,:,:);
        dec_acc_frames_bycl_shuff2 = dec_acc_frames_bycl_shuff(dec_idx,:,:);

        if ~isempty(dec_acc_frames_bycl2)
            figure; hold on; axis tight
            pl_all = cell(num_tn+1);
            for n_tt = 1:num_tn
                %plot(plot_t, dec_acc_frames_bycl_shuff2', color=[0 0 0 0.2])
                pl_all{num_tn+1} = plot(plot_t, mean(dec_acc_frames_bycl_shuff2(:,:,n_tt),1), color=[.5 .5 .5], LineWidth=2);
            end
            for n_tt = 1:num_tn
                col2 = ops.context_types_all_colors2{tn_all(n_tt)};
                %plot(plot_t, dec_acc_frames_bycl2', color=[0 0.4470 0.7410 0.2])
                pl_all{n_tt} = plot(plot_t, mean(dec_acc_frames_bycl2(:,:,n_tt),1), color=col2, LineWidth=2);
            end
            title(sprintf('%s; %s decoder, freqs', reg_labels{n_reg}, params.decoder_type), 'interpreter', 'none');
            legend([pl_all{:}], [ops.context_types_labels(tn_all); {'Shuffle'}]);
            xlim([-0.5, 2.5]);
        end
    end
end


end
