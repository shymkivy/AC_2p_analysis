function [peak_val, peak_loc] = f_get_trial_peak(trial_data_sort, peak_size)

if ~exist('peak_size', 'var')
    peak_size = 1;
end
[num_cells, num_t, ~] = size(trial_data_sort);

peak_pad_left = floor((peak_size - 1)/2);
peak_pad_right = ceil((peak_size - 1)/2);

trial_ave = mean(trial_data_sort,3);
if peak_size > 1
    peak_val = zeros(num_cells, 1);
    peak_loc = zeros(num_cells, 1);
    for n_cell = 1:num_cells
        [~, peak_loc(n_cell)] = max(trial_ave(n_cell,:));
        pk_l = peak_loc(n_cell) - peak_pad_left;
        pk_r = peak_loc(n_cell) + peak_pad_right;
        if pk_l < 1
            pk_l2 = pk_l + 1 - pk_l;
            pk_r2 = pk_r + 1 - pk_l;
        elseif pk_r > num_t
            pk_l2 = pk_l + num_t - pk_r;
            pk_r2 = pk_r + num_t - pk_r;
        else
            pk_l2 = pk_l;
            pk_r2 = pk_r;
        end
        peak_val(n_cell) = mean(trial_ave(n_cell,pk_l2:pk_r2));
    end
else
    [peak_val, peak_loc] = max(trial_ave,[],2);
end


end