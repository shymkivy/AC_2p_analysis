function f_dred_plot_factors(dred_list,trial_types_dred, test_data_ind)



method1 = dred_list(1).method;
num_comp = dred_list(1).n_comp;
num_cv = dred_list(1).cv_num_folds;
colors1 = {'r', 'b', 'g', 'c','m', 'y'};

for n_comp = 1:num_comp
    for n_cv = 1:num_cv
        % first reconstruct for each rank
        d_factors = dred_list([dred_list.n_cv] == n_cv).dred_factors;
        if strcmpi(method1, 'svd')
        elseif strcmpi(method1, 'nmf')
            rank1 = d_factors.d_W(:,n_comp) * d_factors.d_H(n_comp,:);
            rank1_sort = reshape(rank1, size(d_factors.d_W,1),[],sum(~test_data_ind(:,n_cv)));
        elseif strcmpi(method1, 'tca')
            rank1_sort = double(full(ktensor(d_factors.t_factors.lambda(n_comp),...
                                        d_factors.t_factors.U{1}(:,n_comp),...
                                        d_factors.t_factors.U{2}(:,n_comp),...
                                        d_factors.t_factors.U{3}(:,n_comp))));
            %rank1_sort = rank1_sort + d_factors.means;
        end
        trial_types2 = trial_types_dred(~test_data_ind(:,n_cv));
        tt1 = unique(trial_types2);
        %rank1_ave = f_mpl_trial_average(rank1_sort,trial_types2, tt1, 'none');
        
        ymin = 0;
        ymax = 0;
        sp = cell(numel(tt1),1);
        figure;
        for n_tt= 1:numel(tt1)
            sp{n_tt} = subplot(2,3,n_tt); hold on;
            %plot(mean(rank1_sort(:,:,trial_types2 == tt1(n_tt)),3)', 'k', 'LineWidth', 0.1);
            mean_trace = mean(mean(rank1_sort(:,:,trial_types2 == tt1(n_tt)),3),1);
            ymin = min([ymin mean_trace(:)']);
            ymax = max([ymax mean_trace(:)']);
            plot(mean_trace, colors1{n_tt}, 'LineWidth', 2);
            title(num2str(tt1(n_tt)))
        end
        suptitle(sprintf('rank %d, %s, cv %d', n_comp, method1, n_cv))
        for n_tt= 1:numel(tt1)
            ylim(sp{n_tt}, [ymin ymax]);
        end
    end
end

end