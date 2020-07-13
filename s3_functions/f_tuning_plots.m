function f_tuning_plots(data, ops)
    responsive_cells = zeros(numel(ops.conditions_to_analyze),1);
    ctx_tuned_totals = zeros(numel(ops.conditions_to_analyze),10);
    for n_cond = 1:numel(ops.conditions_to_analyze)
        cond_name = ops.conditions{ops.conditions_to_analyze(n_cond)};
        ctx_cells = f_concat_dsets(data.(cond_name).proc, 'context_resp_cells', 1);
        ctx_cells = ctx_cells(:);
        responsive_cells(n_cond) = 100*sum(ctx_cells)/numel(ctx_cells);
        ctx_tuned_totals(n_cond,:) = sum(f_concat_dsets(data.(cond_name).proc, 'tuned_cells_totals', 1),1);

        figure;
        bar(ctx_tuned_totals(n_cond,:), 'FaceColor', [83, 177, 232]/256, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5)
        title([ops.paradigm_type ' ' cond_name 'stimulus tuned cells']);
        
    end
    figure;
    bar(categorical(ops.conditions(ops.conditions_to_analyze)), responsive_cells, 'FaceColor', [83, 177, 232]/256, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5)
    title(sprintf('%s; %s; percent responsive cells, z = %d',ops.paradigm_type, cond_name, ops.z_threshold));
    
    figure;
    bar(sum(ctx_tuned_totals,1), 'FaceColor', [83, 177, 232]/256, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5)
    title([ops.paradigm_type ': All stimulus tuned cells']);
end