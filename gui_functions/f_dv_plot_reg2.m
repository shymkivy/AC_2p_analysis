function f_dv_plot_reg2(app)


%%
region_means_all = cat(3,reg_data.wf_region_means);
region_means_all_in = region_means_all*2^interp_k - 2^interp_k + 1;

num_dset = numel(reg_data);
num_regions = size(reg_data(1).wf_mapping_regions_coords,2);


anchor_means = region_means_all(:,:,anchor_dset);
anchor_im_size = size(reg_data(anchor_dset).wf_mapping_im{freq_im_num})*2^interp_k - 2^interp_k + 1;

%% for current mouse plot wf im
mouse_tag = app.ddata.mouse_tag{1};
current_rdata_idx = strcmpi(mouse_tag,{reg_data.mouse_tag});
current_rdata = reg_data(current_rdata_idx);
current_means = region_means_all(:,:,current_rdata_idx);
current_means_in = region_means_all_in(:,:,current_rdata_idx);
%rdata.regions(1).fov_fname{2}
%rdata.regions(1).region_name{1};

image_wf = current_reg.wf_mapping_im{freq_im_num};
image_wf = image_wf - min(image_wf(:));
image_wf = image_wf/max(image_wf(:));
image_wf_in = interp2(image_wf, interp_k); %

%figure; hold on; imagesc(image_wf_in); 

% register to anchor
current_tform = fitgeotrans(current_means,anchor_means ,'nonreflectivesimilarity');
current_tform.T(3,1:2) = current_tform.T(3,1:2)*2^interp_k - 2^interp_k + 1;
if anchor_dset
    image_wf_in = imwarp(image_wf_in,current_tform ,'OutputView',imref2d(anchor_im_size)); %
    current_means_in2 = current_tform.T' * [current_means_in, ones(size(current_means_in,1),1)]';
    current_means_in3 = current_means_in2(1:2,:)';
else
    current_means_in3 = current_means_in;
end

comb_im_fov = image_wf_in;
comb_im_A = image_wf_in;

if plot_fov
    figure;
    comb_im_fov_fig = imagesc(comb_im_fov); hold on; axis equal tight;
    title(['mouse ' mouse_tag ' fov ' app.ContoursButtonGroup.SelectedObject.Text], 'interpreter', 'none');
    %plot(current_means_in3(:,1),current_means_in3(:,2), 'or')
end

if plot_cell_contours
    figure;
    comb_im_A_fig = imagesc(comb_im_A); hold on; axis equal tight;
    title(['mouse ' mouse_tag ' cell contours ' app.ContoursButtonGroup.SelectedObject.Text], 'interpreter', 'none');
    %plot(current_means_in3(:,1),current_means_in3(:,2), 'or')
end

if app.NewplotsCheckBox.Value
    app.gui_plots.registration_fig = figure;
else
    if isgraphics(app.gui_plots.registration_fig)
        figure(app.gui_plots.registration_fig);
        clf(app.gui_plots.registration_fig);
    else
        app.gui_plots.registration_fig = figure;
    end
end
if ~white_bkg
    imagesc(comb_im_fov); 
end
hold on; axis equal tight;
title(['mouse ' mouse_tag ' tuning ' app.ContoursButtonGroup.SelectedObject.Text], 'interpreter', 'none');

%%

for n_ms = 1:numel(data_mouse_tag)
    mouse_tag2 = data_mouse_tag{n_ms};
    idx_mdata = strcmpi(mouse_tag2,app.data.mouse_tag);
    mdata = app.data(idx_mdata,:);
    
    %%
    current_rdata_idx = strcmpi(mouse_tag2,{reg_data.mouse_tag});
    current_rdata = reg_data(current_rdata_idx);
    rdata = current_rdata.regions;
    
    current_means = region_means_all(:,:,current_rdata_idx);
%     current_means_in = region_means_all_in(:,:,current_rdata_idx);
    
    current_tform_wf = fitgeotrans(current_means,anchor_means ,'nonreflectivesimilarity');
    current_tform_wf.T(3,1:2) = current_tform_wf.T(3,1:2)*2^interp_k - 2^interp_k + 1;
%     if anchor_dset
%         image_wf_in = imwarp(image_wf_in,current_tform_wf ,'OutputView',imref2d(anchor_im_size)); %
%         current_means_in2 = current_tform_wf.T' * [current_means_in, ones(size(current_means_in,1),1)]';
%         current_means_in3 = current_means_in2(1:2,:)';
%     else
%         current_means_in3 = current_means_in;
%     end

%     current_wf = app.data_all(n_ms).wf_mapping_im{n_fr};
% 
%     if app.RegisteronButton.Value
%         current_tform = app.wf_tform_all{n_ms};
%         movingRegistered = imwarp(current_wf,current_tform ,'OutputView',imref2d(size(app.wf_axes_map_mouse.CData))); %
%         %comb_im(movingRegistered>0) = movingRegistered(movingRegistered>0);
%         app.wf_axes_map_mouse.CData = movingRegistered;
%     else
%         app.wf_axes_map_mouse.CData = current_wf;
%     end
% 
    
    for n_dset = 1:size(mdata,1)
        mdata2 = mdata(n_dset,:);
        idx_r = strcmpi(mdata2.area, [rdata.region_name]);
        rdata2 = rdata(:,idx_r);
        stats1 = mdata2.stats{n_pl};
        %%
        if app.ConverttoZCheckBox.Value
            pop_mean_val = stats1.pop_mean_val;
            pop_z_factor = stats1.pop_z_factor;
        end

        if ~isempty(rdata2.regions_tforms)

            %% load images and tform
            image_2p = rdata2.fov_im;

            f = mean(mdata2.OA_data{1}.est.f);
            b = mdata2.OA_data{1}.est.b;
            image_2p_b = reshape(f*b, 256, 256)';

            tform = rdata2.regions_tforms.tform;

            %% load the contours and get coords
            accepted_cells = stats1.accepted_cells;
            A = mdata2.OA_data{1}.est.A(:,accepted_cells);
            A_sum = reshape(mean(A,2), 256, 256)';
            image_A = full(A_sum/max(A_sum(:)));

            num_cells = stats1.num_cells;
            coords = ones(num_cells,3);
            [~, ind1] = max(A);
            [coords(:,1), coords(:,2)] = ind2sub([256 256], ind1);

            %figure; imagesc(image_2p_b)
            %figure; imagesc(A_sum)

            %% test plot
    %         z_thresh = 3;%ops.stat.z_scores_thresh;
    % 
    %         figure;
    %         imagesc(image_A); hold on;
    %         for n_cell = 1:num_cells
    %             [temp_val, ind] = max(onset_fr_peak_mag_ave_z(n_cell,:));
    %             if temp_val > z_thresh
    %                 plot(coords(n_cell,1), coords(n_cell,2), 'o', 'color', ops.context_types_all_colors2{ind},'LineWidth',2);
    %             end
    %             [temp_val, ind] = max(offset_fr_peak_mag_ave_z(n_cell,:));
    %             if temp_val > z_thresh
    %                 plot(coords(n_cell,1), coords(n_cell,2), 'x', 'color', ops.context_types_all_colors2{ind},'LineWidth',2);
    %             end
    %         end
    %         title([reg1.region_name{1} ' ; z=' num2str(z_thresh)]);

            %% interpolate

            image_2p_in = interp2(image_2p, interp_k);
            image_A_in = interp2(image_A, interp_k);

            tform_in = tform;
            tform_in.T(3,1:2) = tform_in.T(3,1:2)*2^interp_k - 2^interp_k + 1;
            
            if anchor_dset
                tform_in.T = tform_in.T * current_tform_wf.T;
            end
            
            coords_in = coords;
            coords_in(:,1:2) = coords_in(:,1:2)*2^interp_k - 2^interp_k + 1;

            %% register 
            if plot_fov
                movingRegistered_fov = imwarp(image_2p_in,tform_in,'OutputView',imref2d(size(image_wf_in)));
                comb_im_fov(movingRegistered_fov>0) = movingRegistered_fov(movingRegistered_fov>0);
                comb_im_fov_fig.CData = comb_im_fov;
            end

            if plot_cell_contours
                movingRegistered_A = imwarp(image_A_in,tform_in,'OutputView',imref2d(size(image_wf_in)));
                comb_im_A(movingRegistered_A>0) = movingRegistered_A(movingRegistered_A>0);
                comb_im_A_fig.CData = comb_im_A;
            end

            coords_tf = (coords_in*tform_in.T)';

            %% plot tuning


            if strcmp(app.ContoursButtonGroup.SelectedObject.Text,'None')
                contour_vals = zeros(num_cells,1);
                use_mag_color_map = 0;
            elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'Tuning type')
                tn_all = f_dv_get_trial_number(app);
                tuning_freq = stats1.peak_val_all(:,tn_all);
                resp_cells = stats1.cell_is_resp(:,tn_all);
                tuning_freq(~resp_cells) = 0;
                [max_val, max_idx] = max(tuning_freq, [], 2);
                contour_vals = max_val;
                contour_mag = max_idx;
                c_lim = [min(max_idx) max(max_idx)];
                use_mag_color_map = 0;
            elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'Tuning mag')
                tn_all = f_dv_get_trial_number(app);
                resp_cells = stats1.cell_is_resp;
                peak_vals = stats1.peak_val_all;
                if app.ConverttoZCheckBox.Value
                    peak_vals = (peak_vals - pop_mean_val)./pop_z_factor;
                end
                peak_vals(~resp_cells) = 0;
                peak_vals2 = max(peak_vals(:,tn_all),[],2);
                contour_vals = peak_vals2;
                contour_mag = peak_vals2; % factor to resize magnitudes to fir color
                c_lim = [min(contour_mag) max(contour_mag)];
                use_mag_color_map = 1;
            elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'SNR')
                SNR_list = mdata2.OA_data{n_pl}.proc.SNR2_vals(accepted_cells);
                contour_mag = SNR_list; % factor to resize magnitudes to fir color
                contour_vals = ones(num_cells,1);
                use_mag_color_map = 1;
            elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'Locomotion')
                resp_cells = stats1.loco_cell;
                if app.ConverttoZCheckBox.Value
                    peak_vals = stats1.loco_z;
                else
                    peak_vals = stats1.loco_corr;
                end
                peak_vals(~resp_cells) = 0;
                contour_vals = peak_vals;
                contour_mag = peak_vals; % factor to resize magnitudes to fir color
                c_lim = [min(contour_mag) max(contour_mag)];
                use_mag_color_map = 1;
            end    

            if use_mag_color_map
                % crop if goes beyond lim
                % create color map
                contour_resolution = 50;
                color_map = jet(round(contour_resolution*c_lim(2)) - round(contour_resolution*c_lim(1))+1); 
                color_index = round(contour_resolution*c_lim(1)):round(contour_resolution*c_lim(2));
            end

            contour_vals(contour_vals<plot_resp_thresh) = 0;

            figure(app.gui_plots.registration_fig);
            for n_cell = 1:num_cells
                if app.SizemagadjustCheckBox.Value
                    if app.ConverttoZCheckBox.Value
                        marker_size2 = marker_size1*contour_vals(n_cell)*app.SizefactorEditField.Value+10;
                    else
                        marker_size2 = marker_size1*contour_vals(n_cell)*app.SizefactorEditField.Value+10;
                    end
                else
                    marker_size2 = marker_size1;
                end
                if contour_vals(n_cell)
                    if use_mag_color_map
                        color1 = color_map(round(contour_resolution*contour_mag(n_cell)) == color_index,:);
                    else
                        color1 = app.ops.context_types_all_colors2{tn_all(contour_mag(n_cell))};
                    end
                    plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', color1,'MarkerSize',marker_size2);
                else
                    if plot_nontuned
                        color1 = [.6 .6 .6];
                        plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', color1,'MarkerSize',marker_size2);
                    end
                end

            end
        end
    end
end
if white_bkg
    app.gui_plots.registration_fig.Children.YDir = 'reverse';
end



end