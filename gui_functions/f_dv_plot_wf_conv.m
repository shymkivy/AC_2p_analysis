function f_dv_plot_wf_conv(app)

n_pl = app.mplSpinner.Value;

interp_k = 0; % number of times splits the coords into two (k)
% if used ad increase fac, 185*2^interp_r_fac - 2^interp_r_fac +1
% num points in between 2^k - 1

marker_size1 = 20;
freq_im_num = app.numfreqSpinner.Value;

anchor_dset = app.anchordsetSpinner.Value;

plot_borders = 0;
plot_nontuned = 0;
white_bkg = 1;
combine_ctx = 1;

%%
tn_all = f_dv_get_trial_number(app);
num_tn = numel(tn_all);

[data1, title_tag] = f_dv_get_data_by_mouse_selection(app);

%% extract all data
num_dsets = size(data1,1);

res_vals_all = cell(num_dsets,1);
resp_cells_all = cell(num_dsets,1);
reg_loc_all = cell(num_dsets,1);
reg_label_all = cell(num_dsets,1);

contour_val = app.ContoursDropDown.Value;

for n_dset = 1:num_dsets
    ddata = data1(n_dset,:);
    
    if ~isempty(ddata.registered_data{n_pl})
        num_cells = ddata.cdata{n_pl}.num_cells;
        stats1 = ddata.stats{n_pl};
        st_mean_mean = stats1.stat_trials_mean_mean;
        st_mean_sem = stats1.stat_trials_mean_sem;
        
        [resp_cells, ~, resp_vals_full] = f_dv_get_resp_vals_cells(app, stats1, tn_all);
        
        res_vals_all{n_dset} = resp_vals_full;
        resp_cells_all{n_dset} = resp_cells;
        reg_loc_all{n_dset} = ddata.registered_data{n_pl}.coords;
        reg_label_all{n_dset} = ddata.registered_data{n_pl}.reg_labels;
    end
end
%%

res_vals_all = cat(1, res_vals_all{:});
resp_cells_all = cat(1, resp_cells_all{:});
reg_loc_all = cat(1, reg_loc_all{:});
reg_label_all = cat(1, reg_label_all{:});

num_cells = size(res_vals_all,1);

res_vals_alln = res_vals_all - min(res_vals_all(:));
res_vals_alln = res_vals_alln./max(res_vals_alln(:));
%%

if strcmpi(contour_val, 'None') || strcmpi(contour_val, 'Tuning type')
    use_mag_color_map = 0;
else
    use_mag_color_map = 1;
end

if use_mag_color_map
    c_lim = [min(resp_cells_all) max(resp_cells_all)];
    % crop if goes beyond lim
    % create color map
    contour_resolution = 10;
    color_map = jet(round(contour_resolution*c_lim(2)) - round(contour_resolution*c_lim(1))+1); 
    color_index = round(contour_resolution*c_lim(1)):round(contour_resolution*c_lim(2));
end

%%
reg_data = app.reg_data;

all_mice = unique(app.data.mouse_tag, 'stable');

anchor_mouse_tag = all_mice{anchor_dset};
anch_idx = strcmpi({reg_data.mouse_tag}, anchor_mouse_tag);
anchor_reg = reg_data(:,anch_idx);

%%
region_means_all = cat(3,reg_data.wf_region_means);

anchor_means = region_means_all(:,:,anchor_dset);
anchor_im_size = size(anchor_reg.wf_mapping_im{freq_im_num})*2^interp_k - 2^interp_k + 1;
%%
current_mouse_tag = app.ddata.mouse_tag{1};
current_idx = strcmpi({reg_data.mouse_tag}, current_mouse_tag);
current_reg = reg_data(:,current_idx);
current_means = region_means_all(:,:,current_idx);


%% for current mouse plot wf im
image_wf = current_reg.wf_mapping_im{freq_im_num};
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

reg_loc_all_in2 = round(reg_loc_all_in);

xymin = min(reg_loc_all_in2);
xymax = max(reg_loc_all_in2);

reg_loc_all_in3 = reg_loc_all_in2 - xymin + 1;
%xymin2 = xymin - xymin + 1;
xymax2 = xymax - xymin + 1;

im_all = cell(num_tn,1);
im_sm_all = cell(num_tn,1);

combine_trials = [18, 28; 20, 30; 29, 19];

for n_tn = 1:num_tn
    im1 = zeros(fliplr(xymax2));
    
    for n_cell = 1:num_cells
        if resp_cells_all(n_cell, n_tn)
            im1(reg_loc_all_in3(n_cell,2), reg_loc_all_in3(n_cell,1)) = res_vals_alln(n_cell,n_tn);
        end
    end
    if combine_ctx
        idx1 = tn_all(n_tn) == combine_trials;
        if sum(idx1(:))
            combine_trials2 = combine_trials(logical(sum(idx1,2)),:);
            tn2 = combine_trials2(~logical(sum(idx1,1)));
            tn3 = (tn_all == tn2);
            for n_cell = 1:num_cells
                if resp_cells_all(n_cell, tn3)
                    im1(reg_loc_all_in3(n_cell,2), reg_loc_all_in3(n_cell,1)) = res_vals_alln(n_cell,tn3);
                end
            end
        end
    end
    
    im1_sm = f_smooth_nd(im1, [5 5]);
    
    im_all{n_tn} = im1;
    im_sm_all{n_tn} = im1_sm;
    
end

imdata1 = cat(1,im_sm_all{:});
max_im_val = max(imdata1(:));

borders = app.border_coords;

for n_tn = 1:num_tn
    color1 = 1 - app.ops.context_types_all_colors2{tn_all(n_tn)};
    im_sm_n = im_sm_all{n_tn}/max_im_val;
    im3 = ones(xymax2(2), xymax2(1), 3).*reshape(color1, 1,1,3);
    im4 = im3.* repmat(im_sm_n, 1, 1, 3);
    figure; 
    subplot(1,3,1);
    imagesc(im_all{n_tn}); axis equal tight; axis off; hold on
    if plot_borders
        for n_reg = 1:4
            pos_in1 = borders{n_reg}.Position;
            pos_in2 = pos_in1*2^interp_k - 2^interp_k + 1;
            plot(pos_in2(:,1)-xymin(1) + 1, pos_in2(:,2)-xymin(2) + 1, '.','Color', app.ops.cond_colors{n_reg}, 'LineWidth', 2);
        end
    end
    subplot(1,3,2);
    imagesc(im_sm_all{n_tn}); axis equal tight; axis off; hold on
    if plot_borders
        for n_reg = 1:4
            pos_in1 = borders{n_reg}.Position;
            pos_in2 = pos_in1*2^interp_k - 2^interp_k + 1;
            plot(pos_in2(:,1)-xymin(1) + 1, pos_in2(:,2)-xymin(2) + 1, '.','Color', app.ops.cond_colors{n_reg}, 'LineWidth', 2);
        end
    end
    subplot(1,3,3);
    imagesc((1 - im4)); axis equal tight; axis off; hold on
    if plot_borders
        for n_reg = 1:4
            pos_in1 = borders{n_reg}.Position;
            pos_in2 = pos_in1*2^interp_k - 2^interp_k + 1;
            plot(pos_in2(:,1)-xymin(1) + 1, pos_in2(:,2)-xymin(2) + 1, '.','Color', app.ops.cond_colors{n_reg}, 'LineWidth', 2);
        end
    end
    
    sgtitle(sprintf('tn %d; cmb ctx = %d', tn_all(n_tn), combine_ctx))
end


triplets_to_plot = [2, 7, 9; 18, 19, 20; 28, 29, 30];
num_trp = size(triplets_to_plot,1);

for n_trp = 1:num_trp
    idx1 = tn_all == triplets_to_plot(n_trp,:)';
    if sum(idx1(:)) == 3
        im3 = zeros(xymax2(2), xymax2(1), 3);
        for n_plt = 1:3
            im3(:,:,n_plt) = im_sm_all{idx1(n_plt,:)};
        end
        im3n = im3/max(im3(:));
        figure; imagesc(im3n*1.2); axis equal tight; hold on;
        if plot_borders
            for n_reg = 1:4
                pos_in1 = borders{n_reg}.Position;
                pos_in2 = pos_in1*2^interp_k - 2^interp_k + 1;
                plot(pos_in2(:,1)-xymin(1) + 1, pos_in2(:,2)-xymin(2) + 1, '.','Color', app.ops.cond_colors{n_reg}, 'LineWidth', 2);
            end
        end
    end
end


%%
for n_tn = 1:num_tn
    f1 = figure;
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
            plot(pos_in2(:,1), pos_in2(:,2), '.','Color', app.ops.cond_colors{n_reg}, 'LineWidth', 2);
        end
    end
    
    color1 = app.ops.context_types_all_colors2{tn_all(n_tn)};
    for n_cell = 1:num_cells
        if resp_cells_all(n_cell, n_tn)
            marker_size2 = res_vals_alln(n_cell,n_tn)*app.SizefactorEditField.Value*marker_size1+5;
            plot(reg_loc_all_in(n_cell,1), reg_loc_all_in(n_cell,2), '.', 'color', color1,'MarkerSize', marker_size2);
        end
    end
    title(sprintf('tn %d', tn_all(n_tn)))
    f1.Children.YDir = 'reverse';
end


end