function f_dv_plot_wf_reg(app)

% app stuff
ops = app.ops;
params = f_dv_gather_params(app);
reg_data = app.reg_data;
reg_fig = app.gui_plots.registration_fig;
borders = app.border_coords;

[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);

n_pl = params.n_pl;

interp_k = 0; % number of times splits the coords into two (k)
% if used ad increase fac, 185*2^interp_r_fac - 2^interp_r_fac +1
% num points in between 2^k - 1

marker_size1 = 20;

plot_borders = params.plot_borders;
plot_nontuned = params.plot_nontuned;
white_bkg = params.white_bkg;

%%
% if strcmpi(app.SelectdatagroupDropDown.Value, 'dataset')
%     data_mouse_tag = app.ddata.mouse_tag;
%     data1 = app.ddata;
% elseif strcmpi(app.SelectdatagroupDropDown.Value, 'mouse')
%     data_mouse_tag = app.ddata.mouse_tag;
%     dset_idx = strcmpi(data_mouse_tag,app.data.mouse_tag);
%     data1 = app.data(dset_idx,:);
% elseif strcmpi(app.SelectdatagroupDropDown.Value, 'all')
%     data_mouse_tag = unique(app.data.mouse_tag, 'stable');
%     
%     dset_idx = false(numel(app.data.mouse_tag),1);
%     for n_ms = 1:numel(data_mouse_tag) 
%         dset_idx2 = strcmpi(data_mouse_tag{n_ms},app.data.mouse_tag);
%         dset_idx = dset_idx + dset_idx2;
%     end
%     dset_idx = logical(dset_idx);
%     data1 = app.data(dset_idx,:);
% end

%% extract all data
num_dsets = size(data,1);

sig_val_all = cell(num_dsets,1);
sig_mag_all = cell(num_dsets,1);
reg_loc_all = cell(num_dsets,1);
reg_label_all = cell(num_dsets,1);

contour_val = params.contours;

for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    
    if ~isempty(ddata.registered_data{n_pl})
        num_cells = ddata.cdata{n_pl}.num_cells;
        stats1 = ddata.stats{n_pl};
        st_mean_mean = stats1.stat_trials_mean_mean;
        st_mean_sem = stats1.stat_trials_mean_sem;


        if strcmpi(contour_val, 'None')
            contour_vals = zeros(num_cells,1);
            contour_mag = zeros(num_cells,1);
        elseif strcmpi(contour_val, 'Tuning type')
            tn_all = f_dv_get_trial_number(params);
            tuning_freq = stats1.peak_vals(:,tn_all);
            resp_cells = stats1.peak_resp_cells(:,tn_all);
            tuning_freq(~resp_cells) = 0;
            [max_val, max_idx] = max(tuning_freq, [], 2);
            contour_vals = max_val;
            contour_mag = max_idx;
        elseif strcmpi(contour_val, 'Tuning magnitude')
            tn_all = f_dv_get_trial_number(params);
            resp_cells = stats1.peak_resp_cells(:,tn_all);
            peak_vals = stats1.peak_vals(:,tn_all);
            if params.convert_to_z
                peak_vals = (peak_vals - st_mean_mean)./st_mean_sem;
            end
            peak_vals(~resp_cells) = 0;
            peak_vals2 = max(peak_vals,[],2);
            contour_vals = peak_vals2;
            contour_mag = peak_vals2; % factor to resize magnitudes to fir color
        elseif strcmpi(contour_val, 'SNR')
            accepted_cells = ddata.cdata{n_pl}.accepted_cells;
            SNR_list = ddata.OA_data{n_pl}.proc.SNR2_vals(accepted_cells);
            contour_mag = SNR_list; % factor to resize magnitudes to fir color
            contour_vals = ones(num_cells,1);
        elseif strcmpi(contour_val, 'Skewness')
            accepted_cells = ddata.cdata{n_pl}.accepted_cells;
            Skew_list = ddata.OA_data{n_pl}.proc.Skewness(accepted_cells);
            contour_mag = Skew_list; % factor to resize magnitudes to fir color
            contour_vals = ones(num_cells,1);
        elseif strcmpi(contour_val, 'Locomotion')
            resp_cells = stats1.loco_cell;
            if params.convert_to_z
                peak_vals = stats1.loco_z;
            else
                peak_vals = stats1.loco_corr;
            end
            peak_vals(~resp_cells) = 0;
            contour_vals = peak_vals';
            contour_mag = peak_vals'; % factor to resize magnitudes to fir color
        end
        sig_val_all{n_dset} = contour_vals;
        sig_mag_all{n_dset} = contour_mag;
        reg_loc_all{n_dset} = ddata.registered_data{n_pl}.coords;
        reg_label_all{n_dset} = ddata.registered_data{n_pl}.reg_labels;
    end
end
%%

sig_val_all = cat(1, sig_val_all{:});
sig_mag_all = cat(1, sig_mag_all{:});
reg_loc_all = cat(1, reg_loc_all{:});
reg_label_all = cat(1, reg_label_all{:});

%%

if strcmpi(contour_val, 'None') || strcmpi(contour_val, 'Tuning type')
    use_mag_color_map = 0;
else
    use_mag_color_map = 1;
end

if use_mag_color_map
    c_lim = [min(sig_mag_all) max(sig_mag_all)];
    % crop if goes beyond lim
    % create color map
    contour_resolution = 10;
    color_map = jet(round(contour_resolution*c_lim(2)) - round(contour_resolution*c_lim(1))+1); 
    color_index = round(contour_resolution*c_lim(1)):round(contour_resolution*c_lim(2));
end

%%

all_mice = unique(data.mouse_tag, 'stable');

anchor_mouse_tag = all_mice{params.anchor_dset};
anch_idx = strcmpi({reg_data.mouse_tag}, anchor_mouse_tag);
anchor_reg = reg_data(:,anch_idx);

%%
region_means_all = cat(3,reg_data.wf_region_means);

anchor_means = region_means_all(:,:,params.anchor_dset);
anchor_im_size = size(anchor_reg.wf_mapping_im{params.freq_im_num})*2^interp_k - 2^interp_k + 1;
%%
current_mouse_tag = ddata.mouse_tag{1};
current_idx = strcmpi({reg_data.mouse_tag}, current_mouse_tag);
current_reg = reg_data(:,current_idx);
current_means = region_means_all(:,:,current_idx);


%% for current mouse plot wf im
image_wf = current_reg.wf_mapping_im{params.freq_im_num};
image_wf = image_wf - min(image_wf(:));
image_wf = image_wf/max(image_wf(:));
image_wf_in = interp2(image_wf, interp_k); %

%figure; hold on; imagesc(image_wf_in); 

% register to anchor
current_tform = fitgeotrans(current_means,anchor_means ,'nonreflectivesimilarity');
current_tform.T(3,1:2) = current_tform.T(3,1:2)*2^interp_k - 2^interp_k + 1;
image_wf_in = imwarp(image_wf_in,current_tform ,'OutputView',imref2d(anchor_im_size)); %

comb_im_fov = image_wf_in;

%%
reg_loc_all_in = reg_loc_all(:,1:2)*2^interp_k - 2^interp_k + 1;

%%

if params.new_plots
    reg_fig = figure(render='painters');
else
    if isgraphics(reg_fig)
        figure(reg_fig);
        clf(reg_fig);
    else
        reg_fig = figure(render='painters');
    end
end
if ~white_bkg
    imagesc(comb_im_fov); 
end
hold on; axis equal tight;
title([title_tag ' tuning ' contour_val], 'interpreter', 'none');

%f1 = figure; %imagesc(im_wf); 
%axis equal tight; hold on;
if plot_borders
    
    for n_reg = 1:4
        pos_in1 = borders{n_reg}.Position;
        pos_in2 = pos_in1*2^interp_k - 2^interp_k + 1;
        plot(pos_in2(:,1), pos_in2(:,2), '.','Color', ops.cond_colors{n_reg}, 'LineWidth', 2);
    end
end
% for n_reg = 1:4
%     reg_idx = reg_label_all == n_reg;
%     plot(reg_loc_all(reg_idx,1),reg_loc_all(reg_idx,2), '.', 'color',app.ops.cond_colors{n_reg});
% end


if plot_nontuned
    for n_cell = 1:size(sig_val_all,1)
        if params.size_mag_adjust
            marker_size2 = marker_size1*sig_val_all(n_cell)*params.size_factor+10;
        else
            marker_size2 = marker_size1;
        end
        if ~sig_val_all(n_cell)
            color1 = [.6 .6 .6];
            plot(reg_loc_all_in(n_cell,1), reg_loc_all_in(n_cell,2), '.', 'color', color1,'MarkerSize',marker_size2);
        end
    end
end
for n_cell = 1:size(sig_val_all,1)
    if params.size_mag_adjust
        marker_size2 = marker_size1*sig_val_all(n_cell)*params.size_factor+10;
    else
        marker_size2 = marker_size1;
    end
    if sig_val_all(n_cell)
        if use_mag_color_map
            color1 = color_map(round(contour_resolution*sig_mag_all(n_cell)) == color_index,:);
        else
            color1 = ops.context_types_all_colors2{tn_all(sig_mag_all(n_cell))};
        end
        plot(reg_loc_all_in(n_cell,1), reg_loc_all_in(n_cell,2), '.', 'color', color1,'MarkerSize',marker_size2);
    end
end

if white_bkg
    reg_fig.Children.YDir = 'reverse';
end


end