function f_plot_trial_indicator2(raster, trial_list, num_bins, ops)
hold on;
[num_cells, num_row] = size(raster);

trial_types = unique(trial_list);
trial_types(trial_types == 0) = [];

num_tt = numel(trial_types);
num_trials = round(num_row/num_bins);

color_seq_tt = zeros(1,num_trials,3);
for n_tt = 1:num_tt
    trials1 = find(trial_list == trial_types(n_tt));
    num_tr = numel(trials1);
    for n_tr_ind = 1:num_tr
        n_tr = trials1(n_tr_ind);
        color_seq_tt(:,n_tr,:) = ops.colors_list{n_tt};
    end
end
color_seq_tt = repmat(color_seq_tt, num_bins ,1, 1);
color_seq_tt = permute(reshape(color_seq_tt, [],1,3),[2 1 3]);

col_width = ceil(num_cells/50);

imagesc(1:num_row,num_cells+(1:col_width),repmat(color_seq_tt,col_width,1,1));
axis tight
end