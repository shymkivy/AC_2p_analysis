function nr_corr = f_mc_compute_nonrigid_motion(Y, params)

[d1, d2, T] = size(Y);

reg_lambda = params.reg_lambda;
block_size = params.nonrigid_block_size;
block_overlap = params.nonrigid_block_overlap;
target_frame = params.target_frame;

num_split = ceil((d1 - block_size)/(block_size-block_overlap) + 1);

cent_coords_m = round(linspace(1, d1, num_split+2));
cent_coords_m2 = cent_coords_m(2:end-1);

cent_coords_n = round(linspace(1, d2, num_split+2));
cent_coords_n2 = cent_coords_n(2:end-1);

coords_mn = zeros(num_split, num_split, 2);
dsall_mn = cell(num_split, num_split);
tic;
fprintf('Block #/%d', num_split^2);
for mm = 1:num_split
    for nn = 1:num_split
        fprintf('..%d', (mm-1)*num_split + nn);
        coords_mn(mm, nn, :) = [cent_coords_m2(mm) cent_coords_n2(nn)];

        m_low = max(round((coords_mn(mm, nn, 1) - block_size/2)),1);
        m_high = min(round((coords_mn(mm, nn, 1) + block_size/2)),d1);

        n_low = max(round((coords_mn(mm, nn, 2) - block_size/2)),1);
        n_high = min(round((coords_mn(mm, nn, 2) + block_size/2)),d2);

        idx_m = m_low:m_high;
        idx_n = n_low:n_high;

        ref_block = target_frame(idx_m, idx_n);
        Y_cut = Y(idx_m, idx_n, :);

        [dsall_mn{mm, nn}] = f_suite2p_reg_compute_cpu(Y_cut, ref_block, reg_lambda);        

    end
end
fprintf('\nDone; compute duration=%.1fsec \n', toc);

% extract the traces into 3d mat 
dsall_m = zeros(num_split,num_split, T);
dsall_n = zeros(num_split,num_split, T);
for mm = 1:num_split
    for nn = 1:num_split
        dsall_m(mm, nn, :) = dsall_mn{mm, nn}(:,1);
        dsall_n(mm, nn, :) = dsall_mn{mm, nn}(:,2);
    end
end

nr_corr.dsall_m = dsall_m;
nr_corr.dsall_n = dsall_n;
nr_corr.cent_coords_n = cent_coords_n;
nr_corr.cent_coords_m = cent_coords_m;
nr_corr.params = params;

if 0
    figure; hold on;
    plot(dsall_mn{1,1})
    plot(dsall_mn{1,2}-5)
    plot(dsall_mn{1,3}-10)
    plot(dsall_mn{1,4}-15)
    title('Some dsall')

%         n_frame = 100;
%         figure; 
%         subplot(3,3,1);
%         imagesc(dsall_m(:,:,n_frame))
%         subplot(3,3,2);
%         imagesc(dsall_m(:,:,n_frame+1))
%         subplot(3,3,3);
%         imagesc(dsall_m(:,:,n_frame+2))
%         subplot(3,3,4);
%         imagesc(dsall_m_sm(:,:,n_frame))
%         subplot(3,3,5);
%         imagesc(dsall_m_sm(:,:,n_frame+1))
%         subplot(3,3,6);
%         imagesc(dsall_m_sm(:,:,n_frame+2))
%         subplot(3,3,7);
%         imagesc(dsall_m_sm_ip(:,:,n_frame))
%         subplot(3,3,8);
%         imagesc(dsall_m_sm_ip(:,:,n_frame+1))
%         subplot(3,3,9);
%         imagesc(dsall_m_sm_ip(:,:,n_frame+2))
end

end