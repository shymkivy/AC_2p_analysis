function f_dv_plot_tuning_dist(app)

plot_ens = app.plotensamblesCheckBox.Value;
frac = app.plotfractionsCheckBox.Value;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

tn_all = f_dv_get_trial_number(app);
num_tn = numel(tn_all);
num_dsets = numel(data.experiment);
reg_all = app.ops.regions_to_analyze;

[region_num, reg_tag] = f_dv_get_region_sel_val(app);

categories = app.ops.context_types_labels(tn_all);

stim_all = zeros(num_tn,1);
for n_dset = 1:num_dsets
    data1 = data(n_dset,:);
    trial_types = data1.trial_types{1};
    
    [~, trial_types_wctx] =  f_s3_add_ctx_trials([], trial_types, app.ddata.MMN_freq{1}, app.ops);

    stim_all = stim_all + sum(trial_types_wctx == app.ops.context_types_all(tn_all)',1)';
end
figure;
bar(categorical(categories,categories), stim_all);
title([title_tag ' ' reg_tag ' stimuli numbers'], 'Interpreter', 'none');
ylabel('number of stimuli');

%%
resp_cell_all = cell(num_dsets,1);
num_cells_all = zeros(num_dsets,numel(region_num));
for n_dset = 1:num_dsets
    data1 = data(n_dset,:);
    
    if ~plot_ens
        stats1 = cat(1,data1.stats{n_pl});
    else
        stats1 = cat(1,data1.ensemble_tuning_stats{n_pl});
    end
    
    num_cells = sum([stats1.num_cells]);
    
    resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], 'Resp split');
    
%     reg_idx = find(strcmpi(reg_all, data1.area));
%     reg_cell_idx = ones(num_cells,1)*reg_idx;
%     if ~isempty(data1.registered_data{1})
%         if app.UseregdatalabelsCheckBox.Value
%             reg_cell_idx = data1.registered_data{1}.reg_labels;
%         end
%     end
    
    if app.UseregdatalabelsCheckBox.Value
        if ~isempty(data1.registered_data{1})
            reg_cell_labels = data1.registered_data{1}.reg_labels;
        else
            reg_cell_labels = zeros(num_cells,1);
        end
    else
        reg_idx = find(strcmpi(reg_all, data1.area));
        reg_cell_labels = ones(num_cells,1)*reg_idx;
    end

    %reg_cell_idx = logical(sum(reg_cell_labels == region_num,2));
    
    cell_is_resp_reg = zeros(num_cells, numel(tn_all), numel(region_num));
    for n_reg = 1:numel(region_num)
        reg_idx = reg_cell_labels == region_num(n_reg);
        cell_is_resp_reg(reg_idx,:,n_reg) = resp_cells(reg_idx,:);
        num_cells_all(n_dset, n_reg) = sum(reg_idx);
    end
    resp_cell_all{n_dset} = cell_is_resp_reg;  
end
resp_cell_all2 = cat(1,resp_cell_all{:});


figure;
if app.poolregionsCheckBox.Value
    resp1= logical(sum(resp_cell_all2,3));
    num_cells1 = sum(num_cells_all,2);
    if frac
        resp2 = sum(resp1,1)/sum(num_cells1);
    else
        resp2 = sum(resp1,1)/num_dsets;
    end
    
    bar(categorical(categories,categories), resp2, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
    title([title_tag ' ' reg_tag  ' tuning distribution; ' num2str(sum(num_cells1)) ' cells'], 'Interpreter', 'none');
    
else
    resp1 = reshape(sum(resp_cell_all2,1), [], numel(reg_all));
    num_cells2 = sum(num_cells_all,1);
    if frac
        resp2 = resp1./num_cells2;
    else
        resp2 = resp1/num_dsets;
    end
    
    br1 = bar(categorical(categories,categories), resp2); %  'EdgeColor',[219, 239, 255]/256, ,'LineWidth',1.5
    title([title_tag ' ' reg_tag  ' tuning distribution'], 'Interpreter', 'none');
    legend(reg_all, 'location', 'northwest');
    for n_br = 1:numel(br1)
        br1(n_br).FaceColor = app.ops.cond_colors{n_br};
    end

end

if frac
    ylabel('Cell fraction');
else
    ylabel('per dset');
end


resp1 = squeeze(logical(sum(resp_cell_all2,2)));

figure; hold on;
if app.poolregionsCheckBox.Value
    num_cells1 = sum(num_cells_all(:));
    resp11 = logical(sum(resp1,2));
    
    if frac
        resp2 = sum(resp11,1)/num_cells1;
    else
        resp2 = sum(resp11,1)/num_dsets;
    end
    
    bar(resp2, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
    title(sprintf('%s %s; fraction of tuned cells; %d cells total', title_tag, reg_tag, sum(num_cells_all(:))), 'Interpreter', 'none');
    
else
    num_cells2 = sum(num_cells_all,1);
    if frac
        resp2 = sum(resp1,1)./num_cells2;
    else
        resp2 = sum(resp1,1)/num_dsets;
    end
    
    bar(categorical(reg_all,reg_all), resp2); %  'EdgeColor',[219, 239, 255]/256, ,'LineWidth',1.5
    title(sprintf('%s %s; fraction of tuned cells; %d cells total', title_tag, reg_tag, sum(num_cells_all(:))), 'Interpreter', 'none');
    for n_br = 1:numel(reg_all)
        br1 = bar(n_br, resp2(n_br));
        br1.FaceColor = app.ops.cond_colors{n_br};
    end
    ylabel('Cell fraction');
    %legend(reg_all, 'location', 'northwest');
end

%%

figure; hold on;
if app.poolregionsCheckBox.Value
    num_cells1 = sum(num_cells_all(:));
    resp11 = logical(sum(resp1,2));
    
    if frac
        resp2 = sum(resp11,1)/num_cells1;
    else
        resp2 = sum(resp11,1)/num_dsets;
    end
    
    bar(resp2, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
    title(sprintf('%s %s; fraction of tuned cells; %d cells total', title_tag, reg_tag, sum(num_cells_all(:))), 'Interpreter', 'none');
    
else
    num_cells2 = sum(num_cells_all,1);
    if frac
        resp2 = sum(resp1,1)./num_cells2;
    else
        resp2 = sum(resp1,1)/num_dsets;
    end
    
    bar(categorical(reg_all,reg_all), resp2); %  'EdgeColor',[219, 239, 255]/256, ,'LineWidth',1.5
    title(sprintf('%s %s; fraction of tuned cells; %d cells total', title_tag, reg_tag, sum(num_cells_all(:))), 'Interpreter', 'none');
    for n_br = 1:numel(reg_all)
        br1 = bar(n_br, resp2(n_br));
        br1.FaceColor = app.ops.cond_colors{n_br};
    end
    ylabel('Cell fraction');
    %legend(reg_all, 'location', 'northwest');
end

%%
if ~app.poolregionsCheckBox.Value
    figure; hold on
    bar(categorical(reg_all,reg_all), num_cells2);
    title(sprintf('%s %s; cell counts; %d cells total', title_tag, reg_tag, sum(num_cells2)), 'Interpreter', 'none');
    ylabel('number of cells');
    for n_br = 1:numel(num_cells2)
        br1 = bar(n_br, num_cells2(n_br));
        br1.FaceColor = app.ops.cond_colors{n_br};
    end
end

%%

loco_cell_reg_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    data1 = data(n_dset,:);
    
    stats1 = cat(1,data1.stats{n_pl});
    num_cells = sum([stats1.num_cells]);
    
    loco_cell = cat(1,[stats1.loco_resp_cells])';
    
    reg_idx = find(strcmpi(reg_all, data1.area));
    reg_cell_idx = ones(num_cells,1)*reg_idx;
    if ~isempty(data1.registered_data{1})
        if app.UseregdatalabelsCheckBox.Value
            reg_cell_idx = data1.registered_data{1}.reg_labels;
        end
    end
    
    loco_cell_reg = zeros(num_cells, numel(region_num));
    for n_reg = 1:numel(region_num)
        reg_idx = reg_cell_idx == region_num(n_reg);
        loco_cell_reg(reg_idx, n_reg) = loco_cell(reg_idx);
    end
    loco_cell_reg_all{n_dset} = loco_cell_reg;  
end

loco_cell2 = cat(1,loco_cell_reg_all{:});

if app.poolregionsCheckBox.Value
    loco_cell3 = sum(loco_cell2,2);
    num_cells1 = sum(num_cells_all,2);
    loco_frac = sum(loco_cell3,1)./sum(num_cells1,1);
    figure; bar(categorical({'All regions'}), loco_frac);
    title([title_tag ' locomotion tuned cells ' reg_tag], 'Interpreter', 'none')
    ylabel('Fraction');
else
    loco_frac = sum(loco_cell2,1)./sum(num_cells_all,1);
    figure; hold on;
    br1 = bar(categorical(reg_all,reg_all), loco_frac);
    for n_br = 1:numel(loco_frac)
        br1 = bar(n_br, loco_frac(n_br));
        br1.FaceColor = app.ops.cond_colors{n_br};
    end
    title([title_tag ' locomotion tuned cells ' reg_tag], 'Interpreter', 'none');
    %legend(reg_all, 'location', 'northwest');
    ylabel('Fraction');

end

%%

%%
%cdata = f_dv_get_cdata(app);
%params = f_dv_gather_params(app);

% %% ensemble stuff
% ens_counts_resp = zeros(num_dsets,1);
% ens_counts_subd = zeros(num_dsets,5);
% ens_totals = zeros(num_dsets,1);
% for n_dset = 1:num_dsets
%     data1 = data(n_dset,:);
%     if ~isempty(data1.ensemble_tuning{n_pl})
%         resp_ens = data1.ensemble_tuning{n_pl}.peak_resp_cells;
%         freq_ens = logical(sum(resp_ens(:,1:10),2));
%         cont_ens = logical(sum(resp_ens(:,[18 28]),2));
%         red_ens = logical(sum(resp_ens(:,[19 29]),2));
%         dev_ens = logical(sum(resp_ens(:,[20 30]),2));
%         loco_ens = logical(data1.ensemble_tuning{1}.loco_cell)';
%         resp_ens2 = logical(freq_ens+ cont_ens+dev_ens+red_ens+loco_ens);
% 
%         ens_counts_resp(n_dset, 1) = sum(resp_ens2);
%         ens_totals(n_dset, 1) = numel(resp_ens2);
% 
%         ens_counts_subd(n_dset, 1) = sum(freq_ens);
%         ens_counts_subd(n_dset, 2) = sum(cont_ens);
%         ens_counts_subd(n_dset, 3) = sum(red_ens);
%         ens_counts_subd(n_dset, 4) = sum(dev_ens);
%         ens_counts_subd(n_dset, 5) = sum(loco_ens);
%     end
% end
% 
% ens_frac = [ens_counts_resp ens_totals-ens_counts_resp]./ens_totals;
% 
% categories = {'stim tuned', 'other'};
% figure;
% bar(categorical(categories,categories), mean(ens_frac,1)); hold on;
% errorbar(categorical(categories,categories), mean(ens_frac,1), std(ens_frac, [], 1)/sqrt(num_dsets-1), '.');
% title('stimulus tuned ensembles')
% ylabel('Fraction')
% 
% ens_frac_subd = [ens_counts_subd ens_totals-ens_counts_resp];
% categories = {'freq', 'cont', 'red', 'dev', 'loco', 'other'};
% figure;
% bar(categorical(categories,categories), mean(ens_frac_subd,1)); hold on;
% errorbar(categorical(categories,categories), mean(ens_frac_subd,1), std(ens_frac_subd, [], 1)/sqrt(num_dsets-1), '.');
% title('stimulus tuned ensembles')
% ylabel('Fraction')
% 
% 
% freq_tuning = zeros(num_dsets,10);
% for n_dset = 1:num_dsets
%     data1 = data(n_dset,:);
%     resp_ens = data1.ensemble_tuning{1}.peak_resp_cells(:,1:10);
%     freq_tuning(n_dset,:) = sum(resp_ens);
% end
% 
% figure;
% bar(mean(freq_tuning,1)); hold on;
% errorbar(mean(freq_tuning,1), std(freq_tuning,[],1)/sqrt(num_dsets-1), '.');
% title('freq tuned ensembles')
% ylabel('Counts')
% 
% cell_counts = zeros(num_dsets,1);
% for n_dset = 1:num_dsets
%     data1 = data(n_dset,:);
%     num_cells = data1.stats{1}.num_cells;
%     cell_counts(n_dset) = num_cells;
% end
% 
% figure; 
% plot(cell_counts, ens_totals, 'o'); hold on;
% if num_dsets >1
%     mdl = fitlm(cell_counts,ens_totals);
%     line_x = min(cell_counts):max(cell_counts);
%     line_y = mdl.Coefficients(1,1).Variables + mdl.Coefficients(2,1).Variables.*line_x;
%     plot(line_x, line_y)
% end
% xlabel('number of cells'); 
% ylabel('number of ensembles');
% title(sprintf('Ensemble num vs cells; %.1f%%', mean(ens_totals./cell_counts)*100))

end
