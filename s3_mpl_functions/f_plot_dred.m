function fig1 = f_plot_dred(dred_data)

sm = 0;

fig1 = {};

method_list = unique({dred_data.method});
num_comp = unique([dred_data.n_comp]);
numFolds = dred_data(1).cv_num_folds;
kernSD = unique([dred_data.kernSD]);
kernSD(kernSD == 0) = [];
method_list_kern = method_list;
method_list_kern(strcmpi(method_list_kern, 'gpfa')) = [];
fin_num = 1;

if numel(kernSD)>1
    fig1{fin_num} = figure;
    fin_num = fin_num + 1;
    subplot(1,2,1); hold on; axis tight;
    plot_n_comp = 5;
    hold on; axis tight;
    [~, ind2] = min(abs(plot_n_comp-num_comp));
    for n_met = 1:numel(method_list_kern)
        ind1 = logical(([dred_data.n_comp] == num_comp(ind2)) .* strcmpi({dred_data.method}, method_list_kern(n_met)));
        
        if sm
            err1 = mean(reshape([dred_data(ind1).train_err_sm],numFolds,[]));
        else
            err1 = mean(reshape([dred_data(ind1).train_err],numFolds,[]));
        end
        
        %plot([dred_data(ind1).sm_kernelSD], [dred_data(ind1).test_err], 'o')
        plot(mean(reshape([dred_data(ind1).kernSD],numFolds,[])),err1, 'LineWidth', 2);
    end
    title(['Training error ' num2str(num_comp(ind2)) ' components']);
    xlabel('smoothing window');
    ylabel('error')
    legend(method_list_kern);
end


if numel(kernSD)>1
    subplot(1,2,2); hold on; axis tight;
    [~, ind2] = min(abs(plot_n_comp-num_comp));
    for n_met = 1:numel(method_list_kern)
        ind1 = logical(([dred_data.n_comp] == num_comp(ind2)) .* strcmpi({dred_data.method}, method_list_kern(n_met)));
        
        if sm
            err1 = mean(reshape([dred_data(ind1).test_err_sm],numFolds,[]));
        else
            err1 = mean(reshape([dred_data(ind1).test_err],numFolds,[]));
        end
        
        %plot([dred_data(ind1).sm_kernelSD], [dred_data(ind1).test_err], 'o')
        plot(mean(reshape([dred_data(ind1).kernSD],numFolds,[])),err1, 'LineWidth', 2);
    end
    title(['Test error ' num2str(num_comp(ind2)) ' components']);
    xlabel('smoothing window');
    ylabel('error')
    legend(method_list_kern);
end


if numel(num_comp)>1
    fig1{fin_num} = figure;
    subplot(1,2,1); hold on; axis tight;
    plot_sm_kernelSD = 200;
    [~, ind2] = min(abs(plot_sm_kernelSD-kernSD));
    for n_met = 1:numel(method_list)
        if sum(strcmpi(method_list(n_met), 'gpfa'))
            kern1 = 0;
        else
            kern1 = kernSD(ind2);
        end
        ind1 = logical(([dred_data.kernSD] == kern1) .* strcmpi({dred_data.method}, method_list(n_met)));
        
        if sm
            err1 = mean(reshape([dred_data(ind1).train_err_sm],numFolds,[]));
        else
            err1 = mean(reshape([dred_data(ind1).train_err],numFolds,[]));
        end
        
        %plot([dred_data(ind1).sm_kernelSD], [dred_data(ind1).test_err], 'o')
        plot(mean(reshape([dred_data(ind1).n_comp],numFolds,[])),err1, 'LineWidth', 2)
    end
    title(['Training error ' num2str(kernSD(ind2)) 'ms kernel SD']);
    xlabel('number of components');
    ylabel('error')
    legend(method_list);
end


if numel(num_comp)>1
    subplot(1,2,2); hold on; axis tight;
    [~, ind2] = min(abs(plot_sm_kernelSD-kernSD));
    for n_met = 1:numel(method_list)
        if strcmpi(method_list(n_met), 'gpfa')
            kern1 = 0;
        else
            kern1 = kernSD(ind2);
        end
        ind1 = logical(([dred_data.kernSD] == kern1) .* strcmpi({dred_data.method}, method_list(n_met)));
        
        if sm
            err1 = mean(reshape([dred_data(ind1).test_err_sm],numFolds,[]));
        else
            err1 = mean(reshape([dred_data(ind1).test_err],numFolds,[]));
        end
        
        %plot([dred_data(ind1).sm_kernelSD], [dred_data(ind1).test_err], 'o')
        plot(mean(reshape([dred_data(ind1).n_comp],numFolds,[])),err1, 'LineWidth', 2)
    end
    title(['Test error ' num2str(kernSD(ind2)) 'ms kernel SD']);
    xlabel('number of components');
    ylabel('error')
    legend(method_list);
end


end