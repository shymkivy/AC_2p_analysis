function [Y_bidi, bidi_out] = f_fix_bidi_shifts3(Y, params)

if ~exist('params', 'var'); params = struct(); end
if ~isfield(params, 'fix_range'); params.fix_range = -20:20; end
if ~isfield(params, 'smooth_std'); params.smooth_std = [0 0.5 2]; end
if ~isfield(params, 'num_iterations'); params.num_iterations = 1; end
if ~isfield(params, 'plot_stuff'); params.plot_stuff = 0; end
if ~isfield(params, 'title_tag'); params.title_tag = 0; end
if ~isfield(params, 'do_interp'); params.do_interp = 1; end
if ~isfield(params, 'normalize'); params.normalize = 1; end

if ~isfield(params, 'laser_open_frac')
    % prairie version 5
    params.laser_open_frac = .65;%.797; % measured 51us up 13us down (says 63 total)28.8; % 14.4
end

smooth_std = params.smooth_std;
fix_range = params.fix_range;
num_iterations = params.num_iterations;
plot_stuff = params.plot_stuff;
title_tag = params.title_tag;
laser_open_frac = params.laser_open_frac;
do_interp = params.do_interp;
normalize = params.normalize;
dens1 = 1;

%%
num_range = numel(fix_range);
[d1, d2, T] = size(Y);

deg_per_fov = 180 * laser_open_frac;
y0 = 1:d1;
z0 = 1:T;

deg0 = linspace(-deg_per_fov/2, deg_per_fov/2, d2/dens1);
x0 = sin(deg0/360*2*pi);
x0n = x0 - min(x0);
x0n = x0n/max(x0n)*(d2-1)+1;


x_coords = 1:dens1:d2;

[X_corr, Y0, Z0] = meshgrid(x_coords, y0, z0);
[X_real, ~, ~] = meshgrid(x0n, y0, z0);   

%%
best_shifts = zeros(T,num_iterations);
best_corr_vals = zeros(T, num_iterations);
best_zero_vals = zeros(T, num_iterations);

idx_1 = 1:2:d1;
idx_2 = 2:2:d1;

Y2 = single(Y);
if do_interp
    % vq = interp1(x,v,xq)
    % Vq = interp2(X,Y,V,Xq,Yq)
    % Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
    fprintf('Interpolating\n');
    Y2 = interp3(X_corr, Y0, Z0, Y2, X_real, Y0, Z0);
end

Y_bidi = Y2;
for n_rep = 1:num_iterations
    fprintf('Bidi fix iter %d; ', n_rep)
    
    Y_odd = Y_bidi(idx_1,:,:);
    Y_even = Y_bidi(idx_2,:,:);
    
    if sum(smooth_std>0)
        tic;
        Y_odd_sm = f_smooth_movie(Y_odd, smooth_std);
        Y_even_sm = f_smooth_movie(Y_even, smooth_std);
        fprintf('smooth duration=%.1fsec; ', toc);
    else
        Y_odd_sm = Y_odd;
        Y_even_sm = Y_even;
    end
    
    tic;
    for n_fr = 1:T
        frame_odd = double(Y_odd_sm(:,:,n_fr));
        
        frame_even = double(Y_even_sm(:,:,n_fr));
        if normalize
            frame_odd = frame_odd/norm(frame_odd, 'fro');
            frame_even = frame_even/norm(frame_even, 'fro');
        end

        corr_vals = zeros(num_range,1);
        max_shift = d2 - max(abs(fix_range));
        for n_sh = 1:num_range
            sh1 = fix_range(n_sh);
            
            start_even = max([1-sh1, 1]);
            end_even = start_even + max_shift-1;
            frame_even_sh = frame_even(:, start_even:end_even);
            start_odd = max([1+sh1, 1]);
            end_odd = start_odd + max_shift-1;
            frame_odd2 = frame_odd(:,start_odd:end_odd);

            if normalize
                frame_even_sh = frame_even_sh/norm(frame_even_sh, 'fro');
                frame_odd2 = frame_odd2/norm(frame_odd2, 'fro');
            end

            corr_val1 = frame_odd2.*frame_even_sh;
            corr_vals(n_sh) = mean(mean(corr_val1));
            if sh1 == 0
                best_zero_vals(n_fr, n_rep) = mean(mean(corr_val1));
            end
        end
        [best_corr_vals(n_fr, n_rep), idx] = max(corr_vals);
        best_shifts(n_fr, n_rep) = fix_range(idx);
%         if 0
%             frame_in = zeros(d1,d2);
%             frame_in(idx_1,:) = frame_odd;
%             frame_in(idx_2,:) = frame_even;            
%             frame_out = f_bidi_apply_shift(frame_in, best_shifts(n_fr, n_rep));
%             raw_in = Y(:,:,n_fr);
%             raw_out = f_bidi_apply_shift(raw_in, best_shifts(n_fr, n_rep));
%             
%             figure; 
%             subplot(2,3,1);imagesc(frame_in); title('frame in');
%             subplot(2,3,2);imagesc(frame_odd); title('odd in');
%             subplot(2,3,3);imagesc(raw_in); title('raw in');
%             subplot(2,3,4);imagesc(frame_out); title('frame out');
%             subplot(2,3,5);imagesc(frame_even); title('even in');
%             subplot(2,3,6);imagesc(raw_out); title('raw out')
%             
%             figure; hold on;
%             plot(fix_range, corr_vals, 'color', [.4 .4 .4]);
%             plot(fix_range(idx), best_corr_vals(n_fr, n_rep), 'or');
%             plot(0, best_zero_vals(n_fr, n_rep), 'og');
%             title(['frame ' num2str(n_fr)])
%             
%         end
        
    end
    fprintf('compute duration=%.1fsec; ', toc);
    
    % apply to the raw
    tic;
    best_shifts_all = sum(best_shifts,2);
    Y_bidi = f_bidi_apply_shift(Y2, best_shifts_all);
    fprintf('apply durration=%.1fsec\n', toc);
    %f_save_mov_YS(Y_bidi, ['C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data\movies\test_bidi_iter' num2str(n_rep) '.h5'], '/mov');
    %f_save_mov_YS(Y_bidi, ['F:\AC_data\test_bidi_iter' num2str(n_rep) '.h5'], '/mov');
    
    
%     k = 14;
%     Y_odd_frame_sh = circshift(Y_odd_frame_ip, round(k), 2);
%     Y_even_frame_sh = circshift(Y_even_frame_ip, -round(k), 2);
%     
%     rgb_im2 = zeros(128,256,3);
%     rgb_im2(:,:,1) = Y_odd_frame_sh;
%     rgb_im2(:,:,2) = Y_even_frame_sh; 
%     rgb_im2 = rgb_im2/max(rgb_im2(:))*1.5;
%     
%     figure; imagesc(rgb_im2)
    
end

% undo  interp
if do_interp
    fprintf('Deinterpolating\n');
    Y_bidi = interp3(X_real, Y0, Z0, Y_bidi, X_corr, Y0, Z0);
end

bidi_out.params = params;
bidi_out.best_shifts = best_shifts;
bidi_out.best_corr_vals = best_corr_vals;
bidi_out.best_zero_vals = best_zero_vals;


if plot_stuff
    
    fr_idx = 2000:3000;
    figure; 
    imagesc(mean(Y_bidi(:,:,fr_idx),3));
    title(sprintf('example ave; laser frac = %.2f', params.laser_open_frac));

    rgb_im = zeros(round(d1/2),d2,3);
    rgb_im(:,:,1) = mean(Y_bidi(idx_1,:,fr_idx),3);
    rgb_im(:,:,2) = mean(Y_bidi(idx_2,:,fr_idx),3);
    rgb_im = rgb_im/max(rgb_im(:))*1.5;
    figure; 
    imagesc(rgb_im)
    title(sprintf('example odd even lines; laser frac = %.2f', params.laser_open_frac));
    
    colors1 = parula(num_iterations);
    figure;
    ax1 = subplot(2,1,1); hold on; axis tight;
    for n_rep = 1:num_iterations
        plot(best_shifts(:,n_rep), 'color', colors1(n_rep,:));
    end
    plot(sum(best_shifts,2), 'k');
    title(['best shift ' title_tag], 'interpreter', 'none');
    ax2 = subplot(2,1,2); hold on; axis tight;
    for n_rep = 1:num_iterations
        plot(best_corr_vals(:,n_rep)-best_zero_vals(:,n_rep), 'color', colors1(n_rep,:));
    end
    %plot(sum(best_corr_vals-best_zero_vals,2), 'k', 'Linewidth', 1);
    title(['best corr -zero ' title_tag], 'interpreter', 'none');
    linkaxes([ax1 ax2], 'x');
end

    
    
    %shift1 = 10;
    %figure; plot(x_coords_ip((1+shift1):end) - x_coords_ip(1:end-shift1))
    
    % first undo the prairie fix
    
    %Y_odd_frame_ip = interp2(X_corr, Y0, Y_odd_frame, X_real, Y0);
    %Y_even_frame_ip = interp2(X_corr, Y0, Y_even_frame, X_real, Y0);
    %Y_target_half_ip = interp2(X_corr, Y0, target_half, X_real, Y0);
    
%     figure;
%     subplot(3,1,1);
%     imagesc(Y_odd_frame_ip);
%     subplot(3,1,2);
%     imagesc(Y_even_frame_ip);
%     subplot(3,1,3);
%     imagesc(Y_target_half_ip);

%
%     target = mean(Y_bidi(:,:,1:500),3);
%     target_half = (target(idx_1,:) + target(idx_2,:))/2;
%     
%     figure; imagesc(target)
%     figure; imagesc(target_half)
%     
%     Y_odd_frame = mean(Y_odd(:,:,2000:3000),3);
%     Y_even_frame = mean(Y_even(:,:,2000:3000),3);
%     
%     x_coord = (1:256) - 256/2;
%     
%     x_coord_tform = x_coord + 0.0005025*x_coord.^2;
%     
%     % 14.4deg
%     
%     figure;
%     s1 = subplot(1,3,1); hold on;
%     imagesc(target_half); axis tight;
%     s1.YDir = 'reverse';
%     title('target frames');
%     s2 = subplot(1,3,2); hold on;
%     imagesc(Y_odd_frame); axis tight;
%     s2.YDir = 'reverse';
%     title('odd row frames');
%     s3 = subplot(1,3,3); hold on;
%     imagesc(Y_even_frame); axis tight;
%     s3.YDir = 'reverse';
%     title('even row frames');
%     
%     pts_all = cell(0,3);
%     pt_idx = 1;
%     num_points = input('how many points to add? (int); 0 to quit:');
%     while num_points
%         for n_pt = 1:num_points
%             rand_pos = [randsample(256, 1), randsample(128, 1)];
%             pts1 = images.roi.Point(s1);
%             pts1.Position = rand_pos;
%             pts1.Color = [1 0 0];
%             pts1.Label = num2str(pt_idx);
%             pts1.LabelAlpha = 0;
% 
%             pts2 = images.roi.Point(s2);
%             pts2.Position = rand_pos;
%             pts2.Color = [1 0 0];
%             pts2.Label = num2str(pt_idx);
%             pts2.LabelAlpha = 0;
%          
%             pts3 = images.roi.Point(s3);
%             pts3.Position = rand_pos;
%             pts3.Color = [1 0 0];
%             pts3.Label = num2str(pt_idx);
%             pts3.LabelAlpha = 0;
%             pt_idx = pt_idx + 1;
%             pts_all = [pts_all; pts1, pts2, pts3];
%         end
%         num_points = input('how many points to add? (int); 0 to quit:');
%     end
%     
%     % collect 
%     
%     %save('bidi_pts.mat', 'pts_all')
%     
%     num_pts = size(pts_all,1);
%     pts_pos_x = zeros(num_pts,3);
%     for n_pt = 1:num_pts
%         pts_pos_x(n_pt,1) = pts_all(n_pt,1).Position(1);
%         pts_pos_x(n_pt,2) = pts_all(n_pt,2).Position(1);
%         pts_pos_x(n_pt,3) = pts_all(n_pt,3).Position(1);
%     end
% 
%         Y_slice = Y_odd_frame(:,130:140);
% 
%     target_halfn = target_half-mean(target_half(:));
%     Y_slicen = Y_slice - mean(Y_slice(:));
%     
%     
%     figure; imagesc(fliplr(Y_slicen))
%     figure; imagesc(target_halfn)
%     
%     x1 = convn(target_halfn, fliplr(Y_slicen));
%     
%     figure; imagesc(x1)
%     
%     figure; plot(sum(x1(126:128,:),1))
%     x_zero = (1:256) - 128;
% 
%     fitobject1 = fit(pts_pos_x(:,1) - 128, pts_pos_x(:,2) - pts_pos_x(:,1),'poly2');
%     fitobject2 = fit(pts_pos_x(:,1) - 128, pts_pos_x(:,3) - pts_pos_x(:,1),'poly2');
%      
%     x_odd_fit = x_zero.^2*fitobject1.p1 + x_zero*fitobject1.p2 + fitobject1.p3;
%     x_even_fit = x_zero.^2*fitobject2.p1 + x_zero*fitobject2.p2 + fitobject2.p3;
%     
%     figure; hold on; axis tight;
%     plot(pts_pos_x(:,1) - 128, pts_pos_x(:,2) - pts_pos_x(:,1), 'o')
%     plot(pts_pos_x(:,1) - 128, pts_pos_x(:,3) - pts_pos_x(:,1), 'o')
%     plot(x_zero, x_odd_fit);
%     plot(x_zero, x_even_fit);
%     xlabel('x odd'); ylabel('x even');
%     title('quad fit');
%     
%     x_odd_fit2 = x_zero - x_odd_fit;
%     x_odd_fit3 = x_odd_fit2 - min(x_odd_fit2);
%     x_odd_fit3 = x_odd_fit3/max(x_odd_fit3) * 255 + 1;
%     
%     x_even_fit2 = x_zero - x_even_fit;
%     x_even_fit3 = x_even_fit2 - min(x_even_fit2);
%     x_even_fit3 = x_even_fit3/max(x_even_fit3) * 255 + 1;
%     
%     x0 = 1:256;
%     y0 = 1:128;
%     
%     [X_odd_fit, ~] = meshgrid(x_odd_fit3, y0);
%     [X_even_fit, ~] = meshgrid(x_even_fit3, y0);
%     [X0, Y0] = meshgrid(x0, y0);
%     
%     Y_odd_frame_ip = interp2(X_odd_fit, Y0, Y_odd_frame, X0, Y0);
%     Y_even_frame_ip = interp2(X_even_fit, Y0, Y_even_frame, X0, Y0);
%     
%     figure;
%     subplot(2,2,1);
%     imagesc(Y_odd_frame);
%     subplot(2,2,2);
%     imagesc(Y_even_frame);
%     subplot(2,2,3);
%     imagesc(Y_odd_frame_ip);
%     subplot(2,2,4);
%     imagesc(Y_even_frame_ip);
%     
%     
%     k = 10;
%     Y_odd_frame_sh = circshift(Y_odd_frame_ip, round(k), 2);
%     Y_even_frame_sh = circshift(Y_even_frame_ip, -round(k), 2);
%     
%     rgb_im = zeros(128,256,3);
%     rgb_im(:,:,1) = Y_odd_frame_sh;
%     rgb_im(:,:,2) = Y_even_frame_sh; 
%     rgb_im = rgb_im/max(rgb_im(:))*1.5;
%     
%     figure; imagesc(rgb_im)
%     
%     im_full = zeros(256, 256+2*k-1);
%     im_full(idx_2,1:256) = Y_even_frame_ip;
%     im_full(idx_1,(1+2*k-1):end) = Y_odd_frame_ip;
%     figure; imagesc(im_full)
%     
    
    
%     prc1 =  99.9;
%     
%     Y_odd_frame2 = Y_odd_frame - min(Y_odd_frame(:));
%     Y_odd_frame2 = Y_odd_frame2/max(Y_odd_frame2(:));
%     val = prctile(Y_odd_frame2(:), prc1);
%     Y_odd_frame2 = Y_odd_frame2/val;
%     Y_odd_frame2(Y_odd_frame2>1) = 1;
%     
%    
%     Y_even_frame2 = Y_even_frame - min(Y_even_frame(:));
%     Y_even_frame2 = Y_even_frame2/max(Y_even_frame2(:));
%     val = prctile(Y_even_frame2(:), prc1);
%     Y_even_frame2 = Y_even_frame2/val;
%     Y_even_frame2(Y_even_frame2>1) = 1;
%     
%     
%     Y_even_frame_sh = circshift(Y_even_frame2, -36, 2);
%     
%     frame_rgb = zeros(size(Y_odd,1), size(Y_odd,2), 3);
%     frame_rgb(:,:,1) = Y_odd_frame2;
%     frame_rgb(:,:,2) = Y_even_frame_sh;
% 
%     figure; imagesc(frame_rgb)
%     
%     colors1 = parula(num_samp);
%     num_samp = 10;
%     pad1 = 20;
%     range1 = (1+pad1):(d2-pad1);
%     
%     Y_odd_frame3 = Y_odd_frame/norm(Y_odd_frame, 'fro');
%     
%     locs1 = round(linspace(1+pad1, d2-pad1, num_samp));
%     cc_list = zeros(numel(range1),1);
%     max_idx_all = zeros(num_samp,1);
%     figure;
%     hold on;
%     for n_samp = 1:num_samp
%         slice1 = Y_even_frame(:,(locs1(n_samp)-pad1):(locs1(n_samp)+pad1));
%         slice2 = slice1/norm(slice1, 'fro');
%         for n_pix = 1:numel(range1)
%             temp = Y_odd_frame3(:,(range1(n_pix)-pad1):(range1(n_pix)+pad1));
%             temp = temp/norm(temp, 'fro');
%             cc_list(n_pix) = sum(sum(temp.*slice2));
%         end
%         plot(cc_list, 'color', colors1(n_samp,:));
%         [~, max_idx_all(n_samp)] = max(cc_list);
%     end
%     
%     
%     Y = max_idx_all(2:end);
%     X = [ones(num_samp-1,1), (2:num_samp)'];
%     b = X\Y;
%     
%     figure; hold on;
%     plot(max_idx_all, 'o-')
%     plot((2:num_samp), X*b)
%     
%     figure; imagesc(slice2)
    

end