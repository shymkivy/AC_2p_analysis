function f_plot_cond_decoding(dec_data_out, params, ops)

cell_range = params.dec_cell_start:params.dec_cell_int:params.dec_cell_max;
max_len = 0;
for n_cond = 1:numel(ops.regions_to_analyze)
    for n_dset = 1:numel(dec_data_out{n_cond})
        max_len = max(max_len, size(dec_data_out{n_cond}{n_dset},1));
    end
end

figure; 
s1= subplot(2,1,1); hold on;
dec_data_means = cell(numel(ops.regions_to_analyze),1);
pl = cell(numel(ops.regions_to_analyze),1);
for n_cond = 1:numel(ops.regions_to_analyze)
    dec_data_means{n_cond} = nan(numel(dec_data_out{n_cond}),max_len);
    for n_dset = 1:numel(dec_data_out{n_cond})
        temp_mean = mean(dec_data_out{n_cond}{n_dset},2);
        dec_data_means{n_cond}(n_dset,1:numel(temp_mean)) = temp_mean;
        if ~isempty(temp_mean)
            p1 = plot(cell_range(1:numel(temp_mean)), temp_mean, 'color', ops.cond_colors{n_cond});
            p1.Color(4) = 0.5;
        end
    end
    pl{n_cond} = plot(cell_range(1:numel(nanmean(dec_data_means{n_cond}))), nanmean(dec_data_means{n_cond}), 'color', ops.cond_colors{n_cond}, 'LineWidth', 2);
end


subplot(2,1,2); hold on;
for n_cond = 1:numel(ops.regions_to_analyze)
    means1 = nanmean(dec_data_means{n_cond});
    sem1 = nanstd(dec_data_means{n_cond})/sqrt(size(dec_data_means{n_cond},1)-1);
    shadedErrorBar_YS(cell_range(1:numel(means1)), means1, sem1, ops.cond_colors{n_cond});
end
legend([pl{:}], ops.regions_to_analyze)
subplot(s1)
end