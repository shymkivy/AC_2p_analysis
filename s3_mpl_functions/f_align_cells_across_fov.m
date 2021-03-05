function f_align_cells_across_fov(data, ops)


mouse_names = unique(data.mouse_tag);

for n_mous = 1:numel(mouse_names)
    dset_idx = strcmpi(mouse_names{n_mous},data.mouse_tag);
    mdata = data(dset_idx,:);
    
    fov_list = unique(mdata.FOV_num);
    
    for n_fov = 1:numel(fov_list)
        
        mdata2 = mdata(mdata.FOV_num == fov_list(n_fov),:);
        
        %%
        num_reg = numel(mdata2.FOV_num);
        
        %% first get the first
        
        n_reg = 1;
        
        f = mean(mdata2.OA_data{n_reg}.est.f);
        b = mdata2.OA_data{n_reg}.est.b;
        im = f*b;
        
        image_2p_b1 = reshape(im,256,256)'/mean(im);%/max(b(:)) /std(im)
        
        figure; imagesc(image_2p_b1); title(['dset ' num2str(n_reg)]);
        
        comp_idx = mdata.OA_data{n_reg}.proc.comp_accepted;
        A1 = mdata.OA_data{n_reg}.est.A(:,comp_idx);
        A1_3d = permute(reshape(full(A1),256,256, []), [2 1 3]);
        A1_proj = sum(A1_3d,3);
        num_cells1 = size(A1,2);
        
        for n_reg = 2:num_reg
            f = mean(mdata2.OA_data{n_reg}.est.f);
            b = mdata2.OA_data{n_reg}.est.b;
            im = f*b;
            image_2p_b2 = reshape(im,256,256)'/mean(im); %/max(b(:)) /std(im)
            
            C = normxcorr2(image_2p_b2,image_2p_b1);
            [ypeak,xpeak] = find(C==max(C(:)));
            yoffSet = ypeak-size(image_2p_b1,1);
            xoffSet = xpeak-size(image_2p_b1,2);

            image_2p_b2_J = imtranslate(image_2p_b2,[xoffSet, yoffSet]);

            %figure; imshowpair(image_2p_b1, image_2p_b2_J,'Scaling','joint');
            
            
            

            comp_idx = mdata.OA_data{n_reg}.proc.comp_accepted;
            A2 = mdata.OA_data{n_reg}.est.A(:,comp_idx);
            A2_3d = permute(reshape(full(A2),256,256, []), [2 1 3]);
            A2_3d = imtranslate(A2_3d,[xoffSet, yoffSet, 0]);
            A2_proj = sum(A2_3d,3);
            num_cells2 = size(A2,2);
            
            c1_c2_corr_err = cell(num_cells1,1);
            
            for n_c1 = 1:num_cells1
                corr_vals = zeros(num_cells2,1);
                for n_c2 = 1:num_cells2
                    corr_vals(n_c2) = mean(mean(A1_3d(:,:,n_c1).*A2_3d(:,:,n_c2)));
                end
                overlap_cells = find(corr_vals>0);
                
                err_vals = zeros(numel(overlap_cells),1);
                for n_c2_idx = 1:numel(overlap_cells)
                    n_c2 = overlap_cells(n_c2_idx);
                    err_vals(n_c2_idx) = mean(mean((A1_3d(:,:,n_c1) - A2_3d(:,:,n_c2)).^2));
                end
                
                c1_c2_corr_err{n_c1} = [ones(numel(overlap_cells),1)*n_c1, overlap_cells, corr_vals(overlap_cells), err_vals];

            end
            
            figure; hold on;
            for n_c = 1:5
                plot(c1_c2_corr_err{n_c}(:,3), c1_c2_corr_err{n_c}(:,4), 'o')
            end
            %figure; imagesc(corr_vals);
            %figure; plot(err_vals);

            figure; imagesc(A1_proj);
            figure; imagesc(A2_proj);
            
            
            f_imagesc_overlay(image_2p_b1, image_2p_b2);
            f_imagesc_overlay(A1_proj, A2_proj);
            
            
        end 
    end
end



end