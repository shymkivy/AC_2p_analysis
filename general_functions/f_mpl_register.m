function [Y_reg, dsall, image_target] = f_mpl_register(Y, params)
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
num_planes = numel(Y);

if isempty(image_target)
    fprintf('Computing to targer image\n');
    image_target = cell(num_planes,1);
    make_image_targe = 1;
else
    make_image_targe = 0;
end
    
dsall = cell(num_planes,num_iterations);
Y_reg = Y;
for n_iter = 1:num_iterations
    for n_pl = 1:num_planes
        fprintf('Registering iter %d plane %d; ', n_iter, n_pl)
        
        %%
        if make_image_targe
            image_target{n_pl} = mean(Y_reg{n_pl},3);
        end
        
        %% smooth movie
        if sum(smooth_std>0)
            tic;
            Y_sm = f_smooth_movie(Y_reg{n_pl}, smooth_std);
            fprintf('smooth duration=%.1fsec; ', toc);
        else
            Y_sm = Y_reg{n_pl};
        end

        %%
        if save_smooth
            temp_fname_sm = sprintf('%s_sm_pl%d_iter%d.h5',save_fname, n_pl, n_iter);
            f_save_mov_YS(Y_sm, [save_path '\' temp_fname_sm], '/mov');
        end
        %%
        tic;
        dsall{n_pl, n_iter} = f_suite2p_reg_compute(Y_sm, image_target{n_pl});
        fprintf('compute duration=%.1fsec; ', toc);
        clear Y_sm;
        
        %%
        tic;
        Y_reg{n_pl} = uint16(f_suite2p_reg_apply(Y_reg{n_pl}, dsall{n_pl}));
        fprintf('apply durration=%.1fsec\n', toc);
        
        if save_reg
            temp_fname_reg = sprintf('%s_reg_pl%d_iter%d.h5',save_fname, n_pl, n_iter);
            f_save_mov_YS(Y_reg{n_pl}, [save_path '\' temp_fname_reg], '/mov');
        end
    end
end
%%
if plot_stuff
    color1 = parula(5);
    figure;
    for n_iter = 1:num_iterations
        subplot(num_iterations*2,1,1 + (n_iter-1)*num_iterations); hold on;
        for n_pl = 1:num_planes
            plot(dsall{n_pl,n_iter}(:,1), 'color', color1(n_pl,:))
        end
        title(sprintf('x shifts, iter %d', n_iter))
        subplot(num_iterations*2,1,2 + (n_iter-1)*num_iterations); hold on;
        for n_pl = 1:num_planes
            plot(dsall{n_pl,n_iter}(:,2), 'color', color1(n_pl,:))
        end
        title(sprintf('y shifts, iter %d', n_iter))
    end
    
end

fprintf('Done\n');

end