function f_mpl_plot_tuning_dset(resp_cells_on, resp_cells_off, dset_params, ops)

%label_depth = 0.003;
context_type_legend = ops.context_types_labels(:);
% plot mmn responsive cells
%ymax = 0;
totals_resp_onset = mean(resp_cells_on);
totals_resp_offset = mean(resp_cells_off);
totals_resp_intersect = mean(resp_cells_on.*resp_cells_off);

figure;
pl1 = if_plot_bar(context_type_legend, totals_resp_onset, totals_resp_offset,totals_resp_intersect, ops);
%ax = gca;
%ymax = max([ymax ax.YLim(2)]);
title(sprintf('%s, dset%d, tuning full; ctxmmn = [%s]', dset_params.cond_name, dset_params.n_dset, num2str(dset_params.ctx_mmn)));
% for n_ctx = 1:6
%     patch([dset_params.ctx_mmn(n_ctx)-0.5, dset_params.ctx_mmn(n_ctx)+0.5 dset_params.ctx_mmn(n_ctx)+0.5, dset_params.ctx_mmn(n_ctx)-0.5],[-label_depth -label_depth -0.001 -0.001], 'r','FaceColor', ops.context_colors{rem(n_ctx-1,3)+1}, 'EdgeColor', 'none');
% end
legend([pl1{:}], {'Onset', 'Offset', 'Both'});
% ylim([-label_depth ymax])


% totals_resp_onset_ctx = totals_resp_onset(:,dset_params.ctx_mmn);
% totals_resp_offset_ctx = totals_resp_offset(:,dset_params.ctx_mmn);
% totals_resp_intersect_ctx = totals_resp_intersect(:,dset_params.ctx_mmn);

% figure;
% pl1 = if_plot_bar(context_type_legend(dset_params.ctx_mmn), totals_resp_onset_ctx, totals_resp_offset_ctx,totals_resp_intersect_ctx, ops);
% title(sprintf('%s, dset%d, tuning context', dset_params.cond_name, dset_params.n_dset));
% legend([pl1{:}], {'Onset', 'Offset', 'Both'});
% ops.context_colors
end

function pl_h = if_plot_bar(categories, resp_on, resp_off, resp_intersect, ops)

pl_h = cell(3,1);
hold on;
pl_h{1} = bar(categorical(categories,categories), resp_on+resp_off-resp_intersect, 'FaceColor', ops.plot_params.color_on, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
pl_h{2} = bar(categorical(categories,categories), resp_off, 'FaceColor', ops.plot_params.color_off, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
pl_h{3} = bar(categorical(categories,categories), resp_intersect, 'FaceColor', ops.plot_params.color_intersect, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);

end