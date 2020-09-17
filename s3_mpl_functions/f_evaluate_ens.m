function ens_eval_out = f_evaluate_ens(ens_out, gt_trials, gt_cells, params)
plot_stuff = f_get_param(params, 'plot_stuff', 0);

%%
if ~isempty(gt_trials)
    clust_tr = ens_out.trials;
    clust_eval_tr = f_evaluate_ens_result(clust_tr.ens_list, gt_trials, plot_stuff);
    if plot_stuff
        suptitle(sprintf('Ens trial detection, mean acc = %.2f', mean(clust_eval_tr.accuracy)));
    end
    %clust_tr.clust_ident = clust_eval_tr.aligned_seq(clust_tr.clust_ident,2);
    fprintf('Ens clustering trials with %.2f accuracy, %d clusters\n', mean(clust_eval_tr.accuracy), numel(clust_tr.clust_label));
    ens_eval_out.clust_eval_tr = clust_eval_tr;
end

%%
if ~isempty(gt_cells)
    clust_cell = ens_out.cells;
    clust_eval_cell = f_evaluate_ens_result(clust_cell.ens_list, gt_cells, plot_stuff);
    if plot_stuff
        suptitle(sprintf('Ens cell detection, mean acc = %.2f', mean(clust_eval_cell.accuracy)));
    end
    %clust_out_cell.clust_ident = clust_eval_cell.ens_perm_ind(clust_cell.clust_ident);
    fprintf('Ens clustering cells with %.2f accuracy, %d clusters\n', mean(clust_eval_cell.accuracy), numel(clust_cell.clust_label));
    ens_eval_out.clust_eval_cell = clust_eval_cell;
end

%%
% 
% 
% ;
%         if compare_ground_truth
%             f_plot_cell_indicator(im1, ens_list_gt, ops);
%             f_plot_trial_indicator2(im1, trial_list_gt, 1, ops);
%         end
%     end
% 
%     pl3_pc = 1:3;
%     %figure; plot3(U(:,1),U(:,2),U(:,3), 'o')
%     f_plot_comp_scatter(U(:,pl3_pc), ens_list_gt);
%     xlabel('pc 1');
%     ylabel('pc 2');
%     zlabel('pc 3');
%     title('U - cells');
% 
% 
%     %figure; plot3(V(:,1),V(:,2),V(:,3), 'o')
%     f_plot_comp_scatter(V(:,pl3_pc), trial_list_gt);
%     xlabel('pc 1');
%     ylabel('pc 2');
%     zlabel('pc 3');
%     title('V - trials')
% end
% 
% if plot_stuff
%     num_plots = ceil(num_LR_comps/3);
%     for n_plt = 1:num_plots
%         dims1 = rem([1 2 3]+(n_plt-1)*3-1, num_LR_comps)+1;
%         f_plot_comp_scatter(X(:,dims1), trial_list_gt);
%         xlabel(sprintf('pc %d', dims1(1)));
%         ylabel(sprintf('pc %d', dims1(2)));
%         zlabel(sprintf('pc %d', dims1(3)));
%         title(sprintf('%s, dset%d, trial type pcs scores(trials)',params.cond_name, params.n_dset));
%     end
% end
% 
% 
% 
end