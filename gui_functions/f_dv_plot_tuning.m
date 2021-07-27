function f_dv_plot_tuning(app)

n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

tn_all = f_dv_get_trial_number(app);

num_dsets = numel(data.experiment);

if strcmpi(app.regiontoplotDropDown.Value, 'all')
    region_num = 0;
elseif strcmpi(app.regiontoplotDropDown.Value, 'A1')
    region_num = 1;
elseif strcmpi(app.regiontoplotDropDown.Value, 'A2')
    region_num = 2;
elseif strcmpi(app.regiontoplotDropDown.Value, 'AAF')
    region_num = 3;
elseif strcmpi(app.regiontoplotDropDown.Value, 'UF')
    region_num = 4;
end

params = f_dv_gather_params(app);

resp_cell_all = zeros(num_dsets, numel(tn_all));
for n_tn = 1:numel(tn_all)
    for n_dset = 1:num_dsets
        fprintf('dset %d\n', n_dset);
        tn1 = tn_all(n_tn);
        
        data1 = data(n_dset,:);
        
        num_cells = data1.stats{1}.num_cells;
        
        if ~region_num
            reg_cell_idx = ones(num_cells,1);
        else
            if ~isempty(data1.registered_data{1})
                reg_cell_idx = data1.registered_data{1}.reg_labels==region_num;
            else
                reg_cell_idx = zeros(num_cells,1);
            end
        end
        
        resp_cells1 = data1.stats{n_pl}.cell_is_resp(:,tn1).*reg_cell_idx;
        
        resp_cell_all(n_dset, n_tn) = sum(resp_cells1);
        
    end
end

categories = app.ops.context_types_labels(tn_all);

if app.poolgroupsCheckBox.Value
    resp1= sum(resp_cell_all,1);
else
    resp1 = resp_cell_all;
end

figure;
bar(categorical(categories,categories), resp1, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
title([title_tag ' ' app.regiontoplotDropDown.Value], 'Interpreter', 'none')


loco_counts = zeros(num_dsets, 4);
no_loco_counts = zeros(num_dsets, 4);
for n_reg = 1:4
    for n_dset = 1:num_dsets
        data1 = data(n_dset,:);
        if ~isempty(data1.registered_data{1})
            reg_cell = data1.registered_data{1}.reg_labels == n_reg;
            loco_cell = data1.stats{1}.loco_cell';
            reg_loco = loco_cell(reg_cell);
            if ~isempty(reg_loco)
                loco_counts(n_dset, n_reg) = sum(reg_loco);
                no_loco_counts(n_dset, n_reg) = sum(~reg_loco);
            end
        end
    end
end

loco_frac = sum(loco_counts,1)./(sum(no_loco_counts,1) + sum(loco_counts,1));

categories = {'A1', 'A2', 'AAF', 'DF'};
figure; bar(categorical(categories,categories), loco_frac);
title([title_tag ' locomotion tuned cells ' app.regiontoplotDropDown.Value], 'Interpreter', 'none')
ylabel('Fraction')



%% ensemble stuff
ens_counts_resp = zeros(num_dsets,1);
ens_counts_subd = zeros(num_dsets,5);
ens_totals = zeros(num_dsets,1);
for n_dset = 1:num_dsets
    data1 = data(n_dset,:);
    resp_ens = data1.ensemble_tuning{1}.cell_is_resp;
    freq_ens = logical(sum(resp_ens(:,1:10),2));
    cont_ens = logical(sum(resp_ens(:,[18 28]),2));
    red_ens = logical(sum(resp_ens(:,[19 29]),2));
    dev_ens = logical(sum(resp_ens(:,[20 30]),2));
    loco_ens = logical(data1.ensemble_tuning{1}.loco_cell)';
    resp_ens2 = logical(freq_ens+ cont_ens+dev_ens+red_ens+loco_ens);
    
    ens_counts_resp(n_dset, 1) = sum(resp_ens2);
    ens_totals(n_dset, 1) = numel(resp_ens2);
    
    ens_counts_subd(n_dset, 1) = sum(freq_ens);
    ens_counts_subd(n_dset, 2) = sum(cont_ens);
    ens_counts_subd(n_dset, 3) = sum(red_ens);
    ens_counts_subd(n_dset, 4) = sum(dev_ens);
    ens_counts_subd(n_dset, 5) = sum(loco_ens);
end

ens_frac = [ens_counts_resp ens_totals-ens_counts_resp]./ens_totals;

categories = {'stim tuned', 'other'};
figure;
bar(categorical(categories,categories), mean(ens_frac,1)); hold on;
errorbar(categorical(categories,categories), mean(ens_frac,1), std(ens_frac, [], 1)/sqrt(num_dsets-1), '.');
title('stimulus tuned ensembles')
ylabel('Fraction')

ens_frac_subd = [ens_counts_subd ens_totals-ens_counts_resp];
categories = {'freq', 'cont', 'red', 'dev', 'loco', 'other'};
figure;
bar(categorical(categories,categories), mean(ens_frac_subd,1)); hold on;
errorbar(categorical(categories,categories), mean(ens_frac_subd,1), std(ens_frac_subd, [], 1)/sqrt(num_dsets-1), '.');
title('stimulus tuned ensembles')
ylabel('Fraction')


freq_tuning = zeros(num_dsets,10);
for n_dset = 1:num_dsets
    data1 = data(n_dset,:);
    resp_ens = data1.ensemble_tuning{1}.cell_is_resp(:,1:10);
    freq_tuning(n_dset,:) = sum(resp_ens);
end

figure;
bar(mean(freq_tuning,1)); hold on;
errorbar(mean(freq_tuning,1), std(freq_tuning,[],1)/sqrt(num_dsets-1), '.');
title('freq tuned ensembles')
ylabel('Counts')

cell_counts = zeros(num_dsets,1);
for n_dset = 1:num_dsets
    data1 = data(n_dset,:);
    num_cells = data1.stats{1}.num_cells;
    cell_counts(n_dset) = num_cells;
end

figure; 
plot(cell_counts, ens_totals, 'o'); hold on;
if num_dsets >1
    mdl = fitlm(cell_counts,ens_totals);
    line_x = min(cell_counts):max(cell_counts);
    line_y = mdl.Coefficients(1,1).Variables + mdl.Coefficients(2,1).Variables.*line_x;
    plot(line_x, line_y)
end
xlabel('number of cells'); 
ylabel('number of ensembles');
title(sprintf('Ensemble num vs cells; %.1f%%', mean(ens_totals./cell_counts)*100))

end
