function f_plot_trial_indicator3(raster, trial_list, num_bins, colors1, xlabels)

hold on;
[num_cells, num_row] = size(raster);

trial_types = unique(trial_list);
trial_types(trial_types == 0) = [];

num_col_row = size(trial_list,2);

num_tt = numel(trial_types);
num_trials = round(num_row/num_bins);

if ~exist('colors1', 'var') || isempty(colors1)
    colors = jet(num_tt);
    for n_tt = 1:num_tt
        colors1{n_tt} = colors(n_tt,:);
    end
end

if ~exist('xlabels', 'var')
    xlabels = [];
end

color_seq_tt = zeros(1, num_col_row,num_trials,3);
for n_tt = 1:num_tt
    for n_r = 1:num_col_row
        trials1 = find(trial_list(:,n_r) == trial_types(n_tt));
        num_tr = numel(trials1);
        for n_tr_ind = 1:num_tr
            n_tr = trials1(n_tr_ind);
            if trial_list(n_tr, n_r)
                color_seq_tt(1, n_r, n_tr, :) = colors1{trial_list(n_tr, n_r)};
            else
                color_seq_tt(1, n_r, n_tr, :) = [0 0 0];
            end
        end
    end
end

col_width = ceil(num_cells/50);

color_seq_tt2 = repmat(color_seq_tt, num_bins*col_width, 1, 1, 1);
% color_seq_tt = permute(reshape(color_seq_tt, num_row, 1, 3), [2 1 3]);
color_seq_tt3 = reshape(color_seq_tt2, num_bins*col_width*num_col_row, num_row, 3);
imagesc(xlabels, num_cells+(1:col_width*num_bins*num_col_row),color_seq_tt3);
axis tight

end