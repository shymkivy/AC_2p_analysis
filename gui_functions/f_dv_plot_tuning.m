function f_dv_plot_tuning(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;
ddata = app.ddata;
mmn_freq = ddata.MMN_freq{1};

%grating_angles = ddata.proc_data{1,1}.stim_params.grating_angles;
grating_angles = [5, 4, 3, 2, 1, 0, -1, -2, -3, -4]*pi/10;

stats1 = app.ddata.stats{n_pl};

tn_all = 1:30;
tn_all_plot = 1:10;
tn_all2 = f_dv_get_trial_number(app);
tn_all_ctx = tn_all2(tn_all2>10);

[resp_cells, ~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], 'Resp split');

bkg_cells_idx = logical(sum(resp_cells(:, tn_all2),2));
bkg_cells_idx(n_cell) = 1;

num_cells = numel(bkg_cells_idx);

if app.ConverttoZCheckBox.Value
    st_mean_mean = stats1.stat_trials_mean_mean;
    st_mean_sem = stats1.stat_trials_mean_sem;
else
    st_mean_mean = zeros(num_cells,1);
    st_mean_sem = ones(num_cells,1);
end

resp_vals2 = (resp_vals - st_mean_mean)./st_mean_sem;

resp_vals_bkg = resp_vals2(bkg_cells_idx,:);
resp_val_cell = resp_vals2(n_cell,:);

z_pad1 = stats1.stat_trials_mean_sem.*stats1.stat_params.z_thresh./st_mean_sem;
base1 = stats1.stat_trials_mean_mean - st_mean_mean;

sem_top_bkg = base1(bkg_cells_idx)+z_pad1(bkg_cells_idx);
sem_bot_bkg = base1(bkg_cells_idx)-z_pad1(bkg_cells_idx);

sem_top_cell = base1(n_cell)+z_pad1(n_cell);
sem_bot_cell = base1(n_cell)-z_pad1(n_cell);

if app.IndividualtrialsCheckBox.Value
    y_lim_max = max([0, max(sem_top_bkg), max(resp_vals_bkg(:))]);
    y_lim_min = min([0, min(sem_bot_bkg), min(resp_vals_bkg(:))]);
else
    y_lim_max = max([0, max(sem_top_cell), max(resp_val_cell(:))]);
    y_lim_min = min([0, min(sem_bot_cell), min(resp_val_cell(:))]);
end
y_lim_max = y_lim_max + y_lim_max*0.05;

y_lim_min = y_lim_min - y_lim_min*0.05;

if app.NewplotsCheckBox.Value
    app.gui_plots.tuning_fig = figure;
else
    if isgraphics(app.gui_plots.tuning_fig)
        figure(app.gui_plots.tuning_fig);
        clf(app.gui_plots.tuning_fig);
    else
        app.gui_plots.tuning_fig = figure;
    end
end


%app.ColoredbytimeCheckBox.Value

%trial_ave_trace = stats1.stat_trials_mean(n_cell,:);
%trial_sem_trace = stats1.stat_trials_sem(n_cell,:);
%stat_window_t = stats1.stat_window_t;
%stat_plot_intsc = logical(logical(sum(stat_window_t'>=plot_t,2)).*logical(sum(stat_window_t'<=plot_t,2)));

cell_is_resp = resp_cells(n_cell,:);
num_tr = numel(tn_all_plot);
if app.PhaseplotCheckBox.Value
    phase1 = grating_angles*2 - pi/2;
    phase2 = repmat(phase1', [1, sum(bkg_cells_idx)]);
    polarplot(0,0);
    hold on;
    if app.IndividualtrialsCheckBox.Value
        if_plot_polar(phase2, resp_vals_bkg(:, tn_all_plot)', [], .5, [.7 .7 .7]);
    end
    if_plot_polar(phase1', resp_val_cell(:, tn_all_plot)', [], 2, [0 0.4470 0.7410]);
    if_plot_polar(phase2, (ones(1, num_tr).*base1(bkg_cells_idx))', 'k', 1);
    if_plot_polar(phase2, (ones(1, num_tr).*sem_top_bkg)', '--k', .5);
    if_plot_polar(phase2, (ones(1, num_tr).*sem_bot_bkg)', '--k', .5);
    rlim([y_lim_min, y_lim_max]);
    for n_tr = 1:10
        if cell_is_resp(n_tr)
            if_plot_polar(phase1(n_tr), resp_val_cell(n_tr), '*g', 2)
        end
    end
    
    if ~isempty(tn_all_ctx)
        for tr_idx = 1:numel(tn_all_ctx)
            n_tr = tn_all_ctx(tr_idx);
            if n_tr <= 20
                x_val = mmn_freq(2);
            elseif n_tr > 20
                x_val = mmn_freq(1);
            end
            if_plot_polar(phase1(x_val), resp_val_cell(:, n_tr)', 'o', 2, app.ops.context_types_all_colors2{n_tr})
        end
    end
else
    hold on; axis tight;
    if app.IndividualtrialsCheckBox.Value
        plot(resp_vals_bkg(:, tn_all_plot)', 'linewidth', .5, 'color', [.7 .7 .7]);
    end
    plot(resp_val_cell(:, tn_all_plot)', 'k', 'linewidth', 2, 'color', [0 0.4470 0.7410])
    plot((ones(1, num_tr).*base1(bkg_cells_idx))', 'k', 'linewidth', 1)
    plot((ones(1, num_tr).*sem_top_bkg)', '--k', 'linewidth', .5)
    plot((ones(1, num_tr).*sem_bot_bkg)', '--k', 'linewidth', .5)
    for n_tr = 1:10
        if cell_is_resp(n_tr)
            plot(n_tr, resp_val_cell(n_tr), '*g', 'linewidth', 2)
        end
    end
    
    if ~isempty(tn_all_ctx)
        for tr_idx = 1:numel(tn_all_ctx)
            n_tr = tn_all_ctx(tr_idx);
            if n_tr <= 20
                x_val = mmn_freq(2);
            elseif n_tr > 20
                x_val = mmn_freq(1);
            end
            plot(x_val, resp_val_cell(:, n_tr)', 'o', 'linewidth', 2, 'color', app.ops.context_types_all_colors2{n_tr})
        end
    end
    
    ylim([y_lim_min, y_lim_max]);
    if app.ConverttoZCheckBox.Value
        ylabel('Z scores');
    else
        ylabel('Normalized response');
    end
end
title(sprintf('Dset %s; Cell %d', app.ddata.dset_name_full{1}, n_cell), 'Interpreter', 'none');

end

%%
function if_plot_polar(phase1, mag1, plotspec1, linewidth1, color1)
    
if exist('plotspec1', 'var') && ~isempty(plotspec1)
    pl1 = polarplot([phase1; phase1(1,:)], [mag1; mag1(1,:)], plotspec1);
else
    pl1 = polarplot([phase1; phase1(1,:)], [mag1; mag1(1,:)]);
end
if exist('linewidth1', 'var') && ~isempty(linewidth1)
    for n_pl = 1:numel(pl1)
    	pl1(n_pl).LineWidth = linewidth1;
    end
end
if exist('color1', 'var') && ~isempty(color1)
    for n_pl = 1:numel(pl1)
    	pl1(n_pl).Color = color1;
    end
end

end

