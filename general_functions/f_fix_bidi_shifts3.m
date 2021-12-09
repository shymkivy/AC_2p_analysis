function [Y_bidi, bidi_out] = f_fix_bidi_shifts3(Y, params)

if ~exist('params', 'var')
    params = struct();
end

if isfield(params, 'fix_range')
    fix_range = params.fix_range;
else
    fix_range = -20:20;
end

if isfield(params, 'smooth_std')
    smooth_std = params.smooth_std;
else
    smooth_std = [0 0.5 2];
end

if isfield(params, 'num_iterations')
    num_iterations = params.num_iterations;
else
    num_iterations = 2;
end

if isfield(params, 'plot_stuff')
    plot_stuff = params.plot_stuff;
else
    plot_stuff = 0;
end

if isfield(params, 'title_tag')
    title_tag = params.title_tag;
else
    title_tag = 0;
end

normalize = 1;

%%
num_range = numel(fix_range);
[d1, d2, T] = size(Y);

best_shifts = zeros(T,num_iterations);
best_corr_vals = zeros(T, num_iterations);
best_zero_vals = zeros(T, num_iterations);

idx_1 = 1:2:d1;
idx_2 = 2:2:d1;

Y_bidi = Y;
for n_rep = 1:num_iterations
    fprintf('Bidi fix iter %d; ', n_rep)
    
    Y_odd = Y_bidi(idx_1,:,:);
    Y_even = Y_bidi(idx_2,:,:);
    
    Y_odd_frame = mean(Y_odd(:,:,2000:end),3);
    Y_even_frame = mean(Y_even(:,:,2000:end),3);
    
    
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
    Y_bidi = f_bidi_apply_shift(Y, best_shifts_all);
    fprintf('apply durration=%.1fsec\n', toc);
    %f_save_mov_YS(Y_bidi, ['C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data\movies\test_bidi_iter' num2str(n_rep) '.h5'], '/mov');
end

bidi_out.params = params;
bidi_out.best_shifts = best_shifts;
bidi_out.best_corr_vals = best_corr_vals;
bidi_out.best_zero_vals = best_zero_vals;

if plot_stuff
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

end