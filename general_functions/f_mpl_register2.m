function [Y_reg, dsall, image_target] = f_mpl_register2(Y, params)
% all inputs and outputs should be 1d cells

if ~exist('params', 'var')
    params = struct();
end

if isfield(params, 'image_target')
    image_target = params.image_target;
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

if isfield(params, 'save_smooth')
    save_smooth = params.save_smooth;
else
    save_smooth = 0;
end

if isfield(params, 'save_reg')
    save_reg = params.save_reg;
else
    save_reg = 0;
end

if isfield(params, 'save_fname')
    save_fname = params.save_fname;
else
    temp_time = clock;
    tag1 = sprintf('%d_%d_%d_%dh_%dm',temp_time(2), temp_time(3), temp_time(1)-2000, temp_time(4), temp_time(5));
    save_fname = ['movie_save_' tag1];
end

if isfield(params, 'save_path')
    save_path = params.save_path;
else
    save_path = '';
end

if isfield(params, 'plot_stuff')
    plot_stuff = params.plot_stuff;
else
    plot_stuff = 0;
end

%%
if isempty(image_target)
    fprintf('Computing target image\n');
    make_image_targe = 1;
else
    make_image_targe = 0;
end

dsall = cell(num_iterations,1);
Y_reg = Y;
for n_iter = 1:num_iterations
    fprintf('Registering iter %d; ', n_iter)

    %%
    if make_image_targe
        image_target = mean(Y_reg, 3);
    end

    %% smooth movie
    if sum(smooth_std>0)
        tic;
        Y_sm = f_smooth_movie(Y_reg, smooth_std);
        fprintf('smooth duration=%.1fsec; ', toc);
    else
        Y_sm = Y_reg;
    end

    %%
    if save_smooth
        temp_fname_sm = sprintf('%s_sm_iter%d.h5',save_fname, n_iter);
        f_save_mov_YS(Y_sm, [save_path '\' temp_fname_sm], '/mov');
    end
    %%
    tic;
    dsall{n_iter} = f_suite2p_reg_compute(Y_sm, image_target);
    fprintf('compute duration=%.1fsec; ', toc);
    clear Y_sm;

    %%
    tic;
    Y_reg = uint16(f_suite2p_reg_apply(Y_reg, dsall{n_iter}));
    fprintf('apply durration=%.1fsec\n', toc);

    if save_reg
        temp_fname_reg = sprintf('%s_reg_iter%d.h5',save_fname, n_iter);
        f_save_mov_YS(Y_reg, [save_path '\' temp_fname_reg], '/mov');
    end
end

%%
if plot_stuff
    figure;
    for n_iter = 1:num_iterations
        subplot(num_iterations,1, n_iter); hold on;
        plot(dsall{n_iter})
        title(sprintf('%s xy shifts, iter %d',save_fname, n_iter), 'interpreter', 'none');
    end
end

end