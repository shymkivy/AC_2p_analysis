function Y_corr = f_mc_apply_nonrigid_corr(Y, corr_data)

[d1, d2, T] = size(Y);

dsall_m = corr_data.dsall_m;
dsall_n = corr_data.dsall_n;
smooth_std = corr_data.params.nonrigid_block_smooth;

cent_coords_m = corr_data.cent_coords_m;
cent_coords_n = corr_data.cent_coords_n;

% process corrections
dsall_m_sm = f_smooth_movie_double(dsall_m, smooth_std);
dsall_m_sm2 = f_pad_matrix(dsall_m_sm, 1, 1);

dsall_n_sm = f_smooth_movie_double(dsall_n, smooth_std);
dsall_n_sm2 = f_pad_matrix(dsall_n_sm, 1, 1);

Y_corr = Y;

block_size = 1000;
num_blocks = ceil(T/block_size);
start1 = 1;

tic;
fprintf('Applying nonrigid reg in blocks; block #/%d:', num_blocks)
for n_bl = 1:num_blocks
    fprintf('..%d', n_bl);
    
    end1 = min((start1 + block_size-1), T);
    block_size2 = (end1 - start1+1);
    
    [X_in, Y_in, Z_in] = meshgrid(single(cent_coords_n), single(cent_coords_m), single(1:block_size2));
    [X_out_bl, Y_out_bl, Z_out_bl] = meshgrid(single(1:d1), single(1:d2), single(1:block_size2));
   
    dsall_m_sm_ip = interp3(X_in, Y_in, Z_in, single(dsall_m_sm2(:,:,start1:end1)), X_out_bl, Y_out_bl, Z_out_bl);
    Y_out2 = Y_out_bl + dsall_m_sm_ip;
    
    dsall_n_sm_ip = interp3(X_in, Y_in, Z_in, single(dsall_n_sm2(:,:,start1:end1)), X_out_bl, Y_out_bl, Z_out_bl);
    X_out2 = X_out_bl + dsall_n_sm_ip;
    
    Y_corr(:,:,start1:end1) = uint16(interp3(single(Y(:,:,start1:end1)), X_out2, Y_out2, Z_out_bl, 'spline', 0));
    
    start1 = end1 + 1;
end
fprintf('\nDone with nonrigid apply; compute time = %.1f\n', toc);

1

% clear dsall_m_sm_ip dsall_n_sm_ip X_out Y_out;
% 
% dsall_m_sm_ip = interp3(X_in, Y_in, Z_in, single(dsall_m_sm2), X_out, Y_out, Z_out);
% Y_out2 = Y_out - dsall_m_sm_ip;
% clear dsall_m_sm_ip;
% 
% dsall_n_sm_ip = interp3(X_in, Y_in, Z_in, single(dsall_n_sm2), X_out, Y_out, Z_out);
% X_out2 = X_out - dsall_n_sm_ip;
% clear dsall_n_sm_ip X_out Y_out;
% 
% Y_corr = uint16(interp3(single(Y), X_out2, Y_out2, Z_out));

end