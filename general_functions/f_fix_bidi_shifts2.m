function [Y_bidi, bidi_out] = f_fix_bidi_shifts2(Y, params)

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
    %figure; hold on;
    for n_fr = 1:T
        frame_lines1 = double(Y_odd_sm(:,:,n_fr));
        frame_lines1 = frame_lines1/norm(frame_lines1, 'fro');
        frame_lines2 = double(Y_even_sm(:,:,n_fr));
        frame_lines2 = frame_lines2/norm(frame_lines1, 'fro');

        corr_vals = zeros(num_range,1);
        for n_sh = 1:num_range
            sh1 = fix_range(n_sh);
            frame_lines2_sh = circshift(frame_lines2, sh1, 2);
            corr_val1 = frame_lines1.*frame_lines2_sh;
            corr_vals(n_sh) = mean(mean(corr_val1));
        end
        [best_corr_vals(n_fr, n_rep), idx] = max(corr_vals);
        best_zero_vals(n_fr, n_rep) = mean(mean(frame_lines1.*frame_lines2));
        best_shifts(n_fr, n_rep) = fix_range(idx);
%         if (best_shifts(n_fr, n_rep))>15
%             frame_in = zeros(d1,d2);
%             frame_in(idx_1,:) = frame_lines1;
%             frame_in(idx_2,:) = frame_lines2;            
%             frame_out = f_bidi_apply_shift(frame_in, best_shifts(n_fr, n_rep));
%             raw_in = Y(:,:,n_fr);
%             raw_out = f_bidi_apply_shift(raw_in, best_shifts(n_fr, n_rep));
%             
%             figure; 
%             subplot(2,3,1);imagesc(frame_in); title('frame in');
%             subplot(2,3,2);imagesc(frame_lines1); title('odd in');
%             subplot(2,3,3);imagesc(raw_in); title('raw in');
%             subplot(2,3,4);imagesc(frame_out); title('frame out');
%             subplot(2,3,5);imagesc(frame_lines2); title('even in');
%             subplot(2,3,6);imagesc(raw_out); title('raw out')
%             
%             figure; hold on;
%             plot(fix_range, corr_vals, 'color', [.4 .4 .4]);
%             plot(fix_range(idx), best_corr_vals(n_fr, n_rep), 'or');
%             plot(0, best_zero_vals(n_fr, n_rep), 'og');
%             title(['frame ' num2str(n_fr)])
%             1;
%         end
%         plot(fix_range, corr_vals, 'color', [.4 .4 .4]);
%         plot(fix_range(idx), best_corr_vals(n_fr, n_rep), 'or');
%         plot(0, best_zero_vals(n_fr, n_rep), 'og');
    end
    fprintf('compute duration=%.1fsec; ', toc);
    
    % apply to the raw
    tic;
    best_shifts_all = sum(best_shifts,2);
    Y_bidi = f_bidi_apply_shift(Y, best_shifts_all);
    fprintf('apply durration=%.1fsec\n', toc);
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