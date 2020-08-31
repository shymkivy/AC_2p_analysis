function f_mpl_plot_tuning_cond(data, ops)

disp('Plotting tuning...')

context_type_legend = ops.context_type_legend(:);
% plot mmn responsive cells
figure;
sp = cell(3,1);
ymax = 0;
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);

    resp_mmn_onset = cat(1,cdata.peak_tuned_trials_onset_ctx{:});
    resp_mmn_offset = cat(1,cdata.peak_tuned_trials_offset_ctx{:});
    totals_resp_mmn_onset = mean(resp_mmn_onset);
    totals_resp_mmn_offset = mean(resp_mmn_offset);
    totals_resp_mmn_intersect = mean(resp_mmn_onset.*resp_mmn_offset);
    
    sp{n_cond} = subplot(ops.plot_params.reg_sm,ops.plot_params.reg_sn,n_cond); 
    if_plot_bar(ops.context_name_full, totals_resp_mmn_onset, totals_resp_mmn_offset,totals_resp_mmn_intersect, ops);
    ax = gca;
    ymax = max([ymax ax.YLim(2)]);
    title(cond_name);
    if n_cond == numel(ops.regions_to_analyze)
        legend('Onset', 'Offset', 'Both');
    end
end
for n_cond = 1:numel(ops.regions_to_analyze)
    ylim(sp{n_cond}, [0 ymax]);
end
suptitle([ops.paradigm_type ' ctx tuned cells']);

% plot freq responsive cells
figure;
ymax = 0;
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    
    resp_freq_onset = cat(1,cdata.peak_tuned_trials_onset{:});
    resp_freq_onset = resp_freq_onset(:,1:10);
    resp_freq_offset = cat(1,cdata.peak_tuned_trials_offset{:});
    resp_freq_offset = resp_freq_offset(:,1:10);
    totals_resp_freq_onset = mean(resp_freq_onset);
    totals_resp_freq_offset = mean(resp_freq_offset);
    totals_resp_freq_intersect = mean(resp_freq_onset.*resp_freq_offset);
    
    sp{n_cond} = subplot(ops.plot_params.reg_sm,ops.plot_params.reg_sn,n_cond); 
    if_plot_bar(context_type_legend(1:10), totals_resp_freq_onset, totals_resp_freq_offset, totals_resp_freq_intersect, ops);
    ax = gca;
    ymax = max([ymax ax.YLim(2)]);
    title(cond_name);
    if n_cond == numel(ops.regions_to_analyze)
        legend('Onset', 'Offset', 'Both');
    end
end
for n_cond = 1:numel(ops.regions_to_analyze)
    ylim(sp{n_cond}, [0 ymax]);
end
suptitle([ops.paradigm_type ' frequency tuned cells']);

% plot freq responsive all freqs
totals_resp_freq_onset = cell(3,1);
totals_resp_freq_offset = cell(3,1);
totals_resp_freq_intersect = cell(3,1);
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    
    resp_freq_onset = cat(1,cdata.peak_tuned_trials_onset{:});
    resp_freq_onset = resp_freq_onset(:,1:10);
    resp_freq_offset = cat(1,cdata.peak_tuned_trials_offset{:});
    resp_freq_offset = resp_freq_offset(:,1:10);
    totals_resp_freq_onset{n_cond} = mean(logical(sum(resp_freq_onset,2)));
    totals_resp_freq_offset{n_cond} = mean(logical(sum(resp_freq_offset,2)));
    totals_resp_freq_intersect{n_cond}= mean(logical(sum(resp_freq_onset,2)).*logical(sum(resp_freq_offset,2)));
end
totals_resp_freq_onset1 = cat(1, totals_resp_freq_onset{:});
totals_resp_freq_offset1 = cat(1, totals_resp_freq_offset{:});
totals_resp_freq_intersect1 = cat(1, totals_resp_freq_intersect{:});
figure;
if_plot_bar(ops.regions_to_analyze, totals_resp_freq_onset1, totals_resp_freq_offset1, totals_resp_freq_intersect1, ops);
title(sprintf('%s; fraction frequency tuned, z = %d',ops.paradigm_type, ops.stat.z_scores_thresh));
legend('Onset', 'Offset', 'Both');


% tunning across regions
% figure;
% bar(categorical(context_type_legend(1:10),context_type_legend(1:10)), sum(freq_onset_tuning_across_reg1,1), 'FaceColor', [83, 177, 232]/256, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5)
% title([ops.paradigm_type ': All freq tuned cells']);


        

end

function if_plot_bar(categories, resp_on, resp_off, resp_intersect, ops)

hold on;
bar(categorical(categories,categories), resp_on+resp_off-resp_intersect, 'FaceColor', ops.plot_params.color_on, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
bar(categorical(categories,categories), resp_off, 'FaceColor', ops.plot_params.color_off, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
bar(categorical(categories,categories), resp_intersect, 'FaceColor', ops.plot_params.color_intersect, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);

end
