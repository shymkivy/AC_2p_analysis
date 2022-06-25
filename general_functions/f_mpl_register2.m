function [Y_reg, dsall, out_frame] = f_mpl_register2(Y, params)
% all inputs and outputs should be 1d cells

if ~exist('params', 'var')
    params = struct();
end

if isfield(params, 'num_iterations')
    num_iterations = params.num_iterations;
else
    num_iterations = 1;
end

if isfield(params, 'smooth_std')
    smooth_std = params.smooth_std;
else
    smooth_std = [0 0 0];
end

if isfield(params, 'reg_lambda')
    reg_lambda = params.reg_lambda;
else
    reg_lambda = [0 0];
end

if isfield(params, 'save_smooth')
    save_all_steps = params.save_all_steps;
else
    save_all_steps = 0;
end

if isfield(params, 'save_fname')
    save_fname = params.save_fname;
else
    temp_time = clock;
    tag1 = sprintf('%d_%d_%d_%dh_%dm',temp_time(2), temp_time(3), temp_time(1)-2000, temp_time(4), temp_time(5));
    save_fname = ['movie_save_' tag1];
end

if isfield(params, 'save_dir')
    save_path = params.save_dir;
else
    save_path = '';
end

if isfield(params, 'plot_stuff')
    plot_stuff = params.plot_stuff;
else
    plot_stuff = 0;
end

%%

num_sm_std = size(smooth_std,1);
num_reg_lambda = size(reg_lambda,1);

dsall = cell(num_iterations,1);
Y_reg = Y;

if save_all_steps
    temp_fname = sprintf('%s_pre_moco.h5',save_fname);
    f_save_mov_YS(Y, [save_path '\' temp_fname], '/mov');
end

for n_iter = 1:num_iterations
    fprintf('Registering iter %d; ', n_iter)

    %%
%     if make_image_targe
%         image_target = mean(Y_reg, 3);
%     end

    %% smooth movie
    
    smooth_std1 = smooth_std(min(n_iter, num_sm_std),:);
    reg_lambda1 = reg_lambda(min(n_iter, num_reg_lambda),:);
   
    if sum(smooth_std1>0)
        tic;
        Y_sm = f_smooth_movie(Y_reg, smooth_std1); % uses ram
        %Y_sm = f_smooth_movie2(Y_reg, smooth_std); % uses GPU
        fprintf('smooth [%.1f %.1f %.1f] duration=%.1fsec; ', smooth_std1(1), smooth_std1(2), smooth_std1(3), toc);
    else
        Y_sm = Y_reg;
    end

    %%
    if save_all_steps
        temp_fname = sprintf('%s_sm_iter%d.h5',save_fname, n_iter);
        f_save_mov_YS(Y_sm, [save_path '\' temp_fname], '/mov');
    end
    %%
    tic;
    [dsall{n_iter}] = f_suite2p_reg_compute(Y_sm, [], reg_lambda1);
    fprintf('compute duration=%.1fsec; ', toc);
    clear Y_sm;
    
    %%
    tic;
    Y_reg = uint16(f_suite2p_reg_apply(Y_reg, dsall{n_iter}));
    fprintf('apply durration=%.1fsec\n', toc);

    if save_all_steps
        temp_fname_reg = sprintf('%s_reg_iter%d.h5',save_fname, n_iter);
        f_save_mov_YS(Y_reg, [save_path '\' temp_fname_reg], '/mov');
    end
end

out_frame = mean(Y_reg,3);

%%
if plot_stuff
    sp_all = cell(num_iterations,1);
    figure;
    for n_iter2 = 1:num_iterations
        sp_all{n_iter2} = subplot(num_iterations,1, n_iter2); hold on;
        plot(dsall{n_iter2})
        title(sprintf('%s xy shifts, iter %d',save_fname, n_iter2), 'interpreter', 'none');
    end
    linkaxes([sp_all{:}], 'x'); axis tight;
end

end