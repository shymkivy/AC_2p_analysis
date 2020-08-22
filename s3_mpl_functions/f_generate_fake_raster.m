function fake_raster_out = f_generate_fake_raster(model_params, source_data)
num_cells = model_params.num_cells;
num_trials = model_params.num_trials;

peak_mag_list = source_data.peak_mag_list;
reliab_list = source_data.reliab_list;
reliab_thresh = model_params.reliab_thresh;
ens_size = model_params.ens_size;
cell_overlap_fraction = model_params.cell_overlap_fraction;
num_corr_trials = model_params.num_corr_trials;

%% make rand raster
raster = zeros(num_cells, num_trials);

% first add random activity
%cell_reliability = linspace(0.15,1,num_cells);
cell_reliability = randsample(reliab_list, num_cells);

active_trials = cell(num_cells,1);
for n_cell = 1:num_cells
    num_events = round(cell_reliability(n_cell)*num_trials);
    chosen_events = randsample(num_trials, num_events);
    raster(n_cell,chosen_events) = randsample(peak_mag_list, num_events);
    active_trials{n_cell} = chosen_events;
end


%% choose ensemble cells 
ens_cell_cell = cell(numel(ens_size),1);
cell_list = find(cell_reliability<=reliab_thresh);
for n_ens = 1:(numel(ens_size)-1)
    num_corr_cells = round(ens_size(n_ens)*cell_overlap_fraction);
    for n_ens2 = (n_ens+1):(numel(ens_size))
        chosen_cells = randsample(cell_list, num_corr_cells);
        ens_cell_cell{n_ens} = cat(1, ens_cell_cell{n_ens}, chosen_cells);
        ens_cell_cell{n_ens2} = cat(1, ens_cell_cell{n_ens2}, chosen_cells);
        cell_list(logical(sum(cell_list==chosen_cells',2))) = [];
    end
end
for n_ens = 1:numel(ens_size)
    chosen_cells = randsample(cell_list,ens_size(n_ens) - numel(ens_cell_cell{n_ens}));
    ens_cell_cell{n_ens} = cat(1, ens_cell_cell{n_ens}, chosen_cells);
    cell_list(logical(sum(cell_list==chosen_cells',2))) = [];
end

%% choose ens trials
ens_trials_cell = cell(numel(ens_size),1);
trial_list = 1:num_trials;
for n_ens = 1:numel(num_corr_trials)
    ens_trials_cell{n_ens} = randsample(trial_list, num_corr_trials(n_ens))';
    trial_list(logical(sum(trial_list == ens_trials_cell{n_ens}))) = [];
end

%% make ens patterns raster
ens_pattern = zeros(num_cells, num_trials);
for n_ens = 1:numel(ens_size)
    cells1 = ens_cell_cell{n_ens};
    trials1 = ens_trials_cell{n_ens};
    for n_cell = 1:numel(cells1)
        ens_pattern(cells1(n_cell),trials1) = 1;
    end
end

peak_mag_list_high = peak_mag_list(peak_mag_list>prctile(peak_mag_list, 30));
ens_pattern_rate = zeros(num_cells, num_trials);
ens_pattern_rate(logical(ens_pattern(:))) = randsample(peak_mag_list_high, sum(ens_pattern(:)));


%% make raster
raster_ens = raster;
num_cell_activations = sum(ens_pattern,2);
% first remove some activations from given raster
for n_cell = 1:num_cells
    if num_cell_activations(n_cell)
        active_bins = find(raster_ens(n_cell,:));
        % remove some active bins to make room for ens bins
        rem_bins = randsample(active_bins, min(numel(active_bins),num_cell_activations(n_cell)));
        raster_ens(n_cell,rem_bins) = 0;
    end
end

% clear ens frames
raster_ens(logical(ens_pattern(:))) = 0;
raster_ens = raster_ens + ens_pattern_rate;

fake_raster_out.raster = raster;
fake_raster_out.raster_ens = raster_ens;
fake_raster_out.ens_cell_cell = ens_cell_cell;
fake_raster_out.ens_trials_cell = ens_trials_cell;
fake_raster_out.ens_pattern_rate = ens_pattern_rate;

end