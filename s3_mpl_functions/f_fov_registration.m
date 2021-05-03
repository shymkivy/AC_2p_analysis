function f_fov_registration(data, ops)

%% for each dset plot things

reg_data_path = 'C:\Users\ys2605\Desktop\stuff\register_2p_to_wf\reg_save.mat';

%% maybe interpolation for better viewing quality of scaled down images
interp_k = 2; % number of times splits the coords into two (k)
% if used ad increase fac, 185*2^interp_r_fac - 2^interp_r_fac +1
% num points in between 2^k - 1

plot_fov = 0;
plot_cell_contours = 0;
plot_freq_tuning_on = 1;
plot_freq_tuning_off = 1;
plot_freq_tuning_combined = 0;
plot_mmn_on = 0;
plot_mmn_off = 0;



%%

reg_data = load(reg_data_path);

mouse_names = unique(data.mouse_tag);

areas_all = ops.regions_to_analyze;


%%
for n_mous = 1:numel(mouse_names)
    dset_idx_m = strcmpi(mouse_names{n_mous},data.mouse_tag);
    mdata = data(dset_idx_m,:);
    
    dset_idx_r = strcmpi(mouse_names{n_mous},{reg_data.data_all.mouse_tag});
    rdata = reg_data.data_all(dset_idx_r);
    
    %rdata.regions(1).fov_fname{2}
    
    %rdata.regions(1).region_name{1};
    
    image_wf = rdata.wf_mapping_im{2};
    %image_wf = rdata.wf_im{1};  
    image_wf = image_wf - min(image_wf(:));
    image_wf = image_wf/max(image_wf(:));
    
    image_wf_in = interp2(image_wf, interp_k);
    
    comb_im_fov = image_wf_in;
    comb_im_A = image_wf_in;
    
    if plot_fov
        figure;
        comb_im_fov_fig = imagesc(comb_im_fov); hold on; axis equal tight;
        title(['mouse ' mdata.mouse_tag{1}], 'interpreter', 'none');
    end
    
    if plot_cell_contours
        figure;
        comb_im_A_fig = imagesc(comb_im_A); hold on; axis equal tight;
        title(['mouse ' mdata.mouse_tag{1}], 'interpreter', 'none');
    end
    
    if plot_freq_tuning_on
        fig1 = figure;
        imagesc(comb_im_fov); hold on; axis equal tight;
        title(['mouse ' mdata.mouse_tag{1} ' onset cells'], 'interpreter', 'none');
    end
    
    if plot_freq_tuning_off
        fig2 = figure;
        imagesc(comb_im_fov); hold on; axis equal tight;
        title(['mouse ' mdata.mouse_tag{1} ' offset cells'], 'interpreter', 'none');
    end
    
    if plot_freq_tuning_combined
        fig3 = figure;
        imagesc(comb_im_fov); hold on; axis equal tight;
        title(['mouse ' mdata.mouse_tag{1}  ' onset+offset cells'], 'interpreter', 'none');
    end
    
    if plot_mmn_on
        fig4 = figure;
        imagesc(comb_im_fov); hold on; axis equal tight;
        title(['mouse ' mdata.mouse_tag{1}  ' mmn on cells'], 'interpreter', 'none');
    end
    
    if plot_mmn_off
        fig5 = figure;
        imagesc(comb_im_fov); hold on; axis equal tight;
        title(['mouse ' mdata.mouse_tag{1}  ' mmn off cells'], 'interpreter', 'none');
    end
    %%
    
    for n_reg = 1:numel(rdata.regions)
        reg1 = rdata.regions(n_reg);
        m_idx = strcmpi(mdata.area, reg1.region_name);
        mdata2 = mdata(m_idx,:);
        if ~isempty(mdata2)
            if ~isempty(rdata.regions(n_reg).regions_tforms)

                %% load images and tform
                image_2p = rdata.regions(n_reg).fov_im;

                f = mean(mdata2.OA_data{1}.est.f);
                b = mdata2.OA_data{1}.est.b;
                image_2p_b = reshape(f*b, 256, 256)';

                tform = rdata.regions(n_reg).regions_tforms.tform;

                %% load the contours and get coords
                comp_idx = mdata2.OA_data{1}.proc.comp_accepted;
                A = mdata2.OA_data{1}.est.A(:,comp_idx);
                A_sum = reshape(mean(A,2), 256, 256)';
                image_A = full(A_sum/max(A_sum(:)));

                num_cells = size(A,2);
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

                %% plot
                
                if plot_freq_tuning_on
                    onset_fr_mag_z = mdata2.tuning_all{1}.peak_tuning_onset.fr_peak_mag_ave_z(:,1:10);
                    figure(fig1);
                    for n_cell = 1:num_cells
                        [temp_val_on, ~] = max(onset_fr_mag_z(n_cell,:));
                        if temp_val_on < ops.stat.z_scores_thresh
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', [.6 .6 .6],'MarkerSize',20);
                        end
                    end
                    for n_cell = 1:num_cells
                        [temp_val_on, ind_on] = max(onset_fr_mag_z(n_cell,:));
                        if temp_val_on > ops.stat.z_scores_thresh
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', ops.context_types_all_colors2{ind_on},'MarkerSize',20);
                        end
                    end
                end
                
                if plot_freq_tuning_off
                    offset_fr_mag_z = mdata2.tuning_all{1}.peak_tuning_offset.fr_peak_mag_ave_z(:,1:10);
                    figure(fig2);
                    for n_cell = 1:num_cells
                        [temp_val_off, ~] = max(offset_fr_mag_z(n_cell,:));
                        if temp_val_off < ops.stat.z_scores_thresh
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', [.6 .6 .6],'MarkerSize',20);
                        end
                    end
                    for n_cell = 1:num_cells
                        [temp_val_off, ind_off] = max(offset_fr_mag_z(n_cell,:));
                        if temp_val_off > ops.stat.z_scores_thresh
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', ops.context_types_all_colors2{ind_off},'MarkerSize',20);
                        end
                    end
                end
                
                if plot_freq_tuning_combined
                    onset_fr_mag_z = mdata2.tuning_all{1}.peak_tuning_onset.fr_peak_mag_ave_z(:,1:10);
                    offset_fr_mag_z = mdata2.tuning_all{1}.peak_tuning_offset.fr_peak_mag_ave_z(:,1:10);
                    figure(fig3);
                    for n_cell = 1:num_cells
                        [temp_val_on, ~] = max(onset_fr_mag_z(n_cell,:));
                        [temp_val_off, ~] = max(offset_fr_mag_z(n_cell,:));
                        if and(temp_val_on < ops.stat.z_scores_thresh, temp_val_off < ops.stat.z_scores_thresh)
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), 'o', 'color', [.6 .6 .6],'LineWidth',2);
                        end
                    end
                    for n_cell = 1:num_cells
                        [temp_val_on, ind_on] = max(onset_fr_mag_z(n_cell,:));
                        [temp_val_off, ind_off] = max(offset_fr_mag_z(n_cell,:));
                        if temp_val_on > ops.stat.z_scores_thresh
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), 'o', 'color', ops.context_types_all_colors2{ind_on},'LineWidth',2);
                        end
                        
                        if temp_val_off > ops.stat.z_scores_thresh
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), 'x', 'color', ops.context_types_all_colors2{ind_off},'LineWidth',2);
                        end
                    end
                end
                
                sum(sum(mdata2.peak_tuned_trials_offset_ctx{1},2));
                
                if plot_mmn_on
                    onset_mmn_mag_z = mdata2.tuning_all{1}.peak_tuning_onset.fr_peak_mag_ave_z(:,mdata2.ctx_mmn{1});
                    figure(fig4);
%                     for n_cell = 1:num_cells
%                         [temp_val_on, ~] = max(onset_mmn_mag_z(n_cell,:));
%                         if temp_val_on < ops.stat.z_scores_thresh
%                             plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', [.6 .6 .6],'MarkerSize',20);
%                         end
%                     end
                    for n_cell = 1:num_cells
                        [temp_val_on, ind_on] = max(onset_mmn_mag_z(n_cell,:));
                        if temp_val_on > ops.stat.z_scores_thresh
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', ops.context_types_all_colors2{mdata2.ctx_mmn{1}(ind_on)},'MarkerSize',20);
                        end
                    end
                    
                end
                
                if plot_mmn_off
                    offset_mmn_mag_z = mdata2.tuning_all{1}.peak_tuning_offset.fr_peak_mag_ave_z(:,mdata2.ctx_mmn{1});
                    figure(fig5);
%                     for n_cell = 1:num_cells
%                         [temp_val_off, ~] = max(offset_mmn_mag_z(n_cell,:));
%                         if temp_val_off < ops.stat.z_scores_thresh
%                             plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', [.6 .6 .6],'MarkerSize',20);
%                         end
%                     end
                    for n_cell = 1:num_cells
                        [temp_val_off, ind_on] = max(offset_mmn_mag_z(n_cell,:));
                        if temp_val_off > ops.stat.z_scores_thresh
                            plot(coords_tf(1,n_cell), coords_tf(2,n_cell), '.', 'color', ops.context_types_all_colors2{mdata2.ctx_mmn{1}(ind_on)},'MarkerSize',20);
                        end
                    end
                    
                end
                
            end
        end
    end

    
end

if logical(plot_freq_tuning_combined+plot_freq_tuning_on+plot_freq_tuning_off)
    leg_im = reshape(cat(2,ops.context_types_all_colors2{1:10})',10,1,3);
    figure; imagesc(leg_im); ylabel('low to high freqs')
end


end