function [Y, best_shifts] = f_fix_bidi_shifts(Y, params)

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

if isfield(params, 'num_reps')
    num_reps = params.num_reps;
else
    num_reps = 2;
end

%%
num_range = numel(fix_range);
[d1, d2, T] = size(Y);

best_shifts = zeros(T,num_reps);
best_corr_vals = zeros(T, num_reps);
best_zero_vals = zeros(T, num_reps);

idx_1 = 1:2:d1;
idx_2 = 2:2:d1;

Y_bidi = Y;
for n_rep = 1:num_reps
    
    if sum(smooth_std>0)
        Y_sm = f_smooth_movie(Y_bidi, smooth_std);
    else
        Y_sm = Y_bidi;
    end
    
    %figure; hold on;
    f = waitbar(0,'Computing bidi shifts...');
    for n_fr = 1:T
        temp_fr = Y_sm(:,:,n_fr);
        frame_lines1 = temp_fr(idx_1,:);
        frame_lines2 = temp_fr(idx_2,:);

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
        
%         plot(fix_range, corr_vals, 'color', [.4 .4 .4]);
%         plot(fix_range(idx), best_corr_vals(n_fr, n_rep), 'or');
%         plot(0, best_zero_vals(n_fr, n_rep), 'og');
        waitbar(n_fr/T,f);
    end
    close(f)
    
    % apply to the raw
    best_shifts_all = sum(best_shifts,2);
    for n_fr = 1:T
        frame_in = Y(:,:,n_fr);
        new_frame = frame_in - frame_in;

        rem1 = rem(best_shifts_all(n_fr),2);
        half_shift = (best_shifts_all(n_fr)-rem1)/2;

        % shift odd lines right
        shift1 = -half_shift;
        start_in1 = max([1-shift1 1]);
        end_in1 = min([d2-shift1, d2]);
        start_out1 = max([1+shift1 1]);
        end_out1 = min([d2+shift1, d2]);
        new_frame(idx_1,start_out1:end_out1) = frame_in(idx_1,start_in1:end_in1);


        % shift even lines left
        shift2 = half_shift + rem1;
        start_in2 = max([1-shift2 1]);
        end_in2 = min([d2-shift2, d2]);
        start_out2 = max([1+shift2 1]);
        end_out2 = min([d2+shift2, d2]);
        new_frame(idx_2,start_out2:end_out2) = frame_in(idx_2,start_in2:end_in2);

        Y_bidi(:,:,n_fr) = new_frame;
    end
end


figure;
ax1 = subplot(2,1,1); hold on; axis tight;
plot(best_shifts);
plot(sum(best_shifts,2), 'k');
title('best shift');
ax2 = subplot(2,1,2); hold on; axis tight;
plot(best_corr_vals-best_zero_vals);
%plot(sum(best_corr_vals-best_zero_vals,2), 'k', 'Linewidth', 1);
title('best corr -zero');
linkaxes([ax1 ax2], 'x')


f_save_mov_YS(Y, ['C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data\movies\test_no_bidi.hdf5'], '/mov')

f_save_mov_YS(Y_bidi, ['C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data\movies\test_bidi.hdf5'], '/mov')

subplot(3,1,3);
plot(best_zero_vals);
title('best zero vals');

for n_rep = 1:num_reps
    subplot(num_reps, 1, n_rep)
    plot(best_shifts(:,n_rep)); title(['repetition ' num2str(n_rep)]);
end


end