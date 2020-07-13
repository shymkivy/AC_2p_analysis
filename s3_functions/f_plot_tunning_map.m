function f_plot_tunning_map(data, ops, area1, dataset_num)

cells_ind = data.(area1).num_dataset ==dataset_num;
trial_ave = data.(area1).trial_ave_z(cells_ind,:,:);
%tuned_cells = data.(area1).tuned_cells(cells_ind);

cell_locs = data.(area1).cell_loc(cells_ind);


MMN_freq = data.(area1).MMN_freq(cells_ind,:);
MMN_freq = MMN_freq(1,:);


plot_cell_loc(trial_ave, MMN_freq, cell_locs, ops, [0, 0]);
title(sprintf('%s, dataset %d, tuning', area1, dataset_num));
% plot_cell_loc(trial_ave, MMN_freq, cell_locs, ops, [1, 0]);
% title(sprintf('%s, dataset %d, dev', area1, dataset_num));
% plot_cell_loc(trial_ave, MMN_freq, cell_locs, ops, [0, 1]);
% title(sprintf('%s, dataset %d, red', area1, dataset_num));
% plot_cell_loc(trial_ave, MMN_freq, cell_locs, ops, [1, 1]);
% title(sprintf('%s, dataset %d, all', area1, dataset_num));

% 
% ops.win.onset_window_frames
% ops.time_stim_window(ops.win.onset_window_frames)
% ops.time_stim_window(ops.win.offset_window_frames)
% 
% plot_cells = 1:10;
% plot_trace = trial_ave(:,10:end,:); %ops.time_stim_window>=0
% ylims1 = [0 max(plot_trace(:))];
% 
% plot_time = ops.time_stim_window(10:end);
% for ii =  plot_cells
%     for jj = 1:10
%         subplot(10, 10, 10*(ii-1)+jj)
%         plot(plot_time, plot_trace(ii,:,jj))
%         ylim(ylims1);
%     end
% end




end


function plot_cell_loc(trial_ave, MMN_freq, cell_locs, ops, dev_red)

full_tunning = squeeze(mean(trial_ave(:,ops.win.onset_window_frames,:),2));

[~, tunning_ind] = max(full_tunning(:,1:10), [], 2);
tuned_cells = logical(sum(full_tunning(:,1:10)>3,2));

dev1 = (full_tunning(:,28)>3);
dev2 = (full_tunning(:,19)>3);
red1 = logical(sum(full_tunning(:,20:27)>3,2));
red2 = logical(sum(full_tunning(:,11:18)>3,2));

%figure; imagesc(tunning(tuned_cells,:))


im_loc = 256*ones(256,256,3);
figure; imagesc(im_loc);
axis square
axis off
jet_tunmap = jet(10);
colormap gray
hold on;
mrk_size = 8;
for n_cell = 1:numel(cell_locs)
    temp_cell = cell_locs{n_cell};
    if tuned_cells(n_cell)
        plot(temp_cell(5,1),temp_cell(5,2), 'o', 'MarkerSize',mrk_size, 'Color', jet_tunmap(tunning_ind(n_cell),:));
    else
        plot(temp_cell(5,1),temp_cell(5,2), 'o', 'MarkerSize',mrk_size, 'Color', [.7 .7 .7]);
    end
    if dev_red(1)
        if dev1(n_cell)
            plot(temp_cell(5,1),temp_cell(5,2), '+', 'MarkerSize',mrk_size,'LineWidth',2, 'Color', jet_tunmap(MMN_freq(1),:));
        end

        if dev2(n_cell)
            plot(temp_cell(5,1),temp_cell(5,2), '+', 'MarkerSize',mrk_size,'LineWidth',2, 'Color', jet_tunmap(MMN_freq(2),:));
        end
    end
    if dev_red(2)
        if red1(n_cell)
            plot(temp_cell(5,1),temp_cell(5,2), 'x', 'MarkerSize',mrk_size,'LineWidth',2, 'Color', jet_tunmap(MMN_freq(1),:));
        end
        if red2(n_cell)
            plot(temp_cell(5,1),temp_cell(5,2), 'x', 'MarkerSize',mrk_size,'LineWidth',2, 'Color', jet_tunmap(MMN_freq(2),:));
        end
    end
end
%jet_tunmap = jet(10);
% figure;
% imagesc(reshape(jet_tunmap,10,1,3))
% ax = gca;
% ax.XTick = [];
% title('freq colormap');


end
