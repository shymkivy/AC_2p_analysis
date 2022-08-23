function [Y_reg, mc_out] = f_mc_rigid(Y, params)
% all inputs and outputs should be 1d cells

if ~exist('params', 'var'); params = struct(); end
if ~isfield(params, 'num_iterations'); params.num_iterations = 1; end
if ~isfield(params, 'smooth_std'); params.smooth_std = [0 0 0]; end
if ~isfield(params, 'reg_lambda'); params.reg_lambda = [0 0]; end
if ~isfield(params, 'high_val_cut_thresh'); params.high_val_cut_thresh = 0; end
if ~isfield(params, 'high_val_cut_per_frame'); params.high_val_cut_per_frame = 1; end

if ~isfield(params, 'save_all_steps'); params.save_all_steps = 0; end
if ~isfield(params, 'save_dir'); params.save_dir = ''; end
if ~isfield(params, 'plot_stuff'); params.plot_stuff = 0; end


if ~isfield(params, 'save_fname')
    temp_time = clock;
    tag1 = sprintf('%d_%d_%d_%dh_%dm',temp_time(2), temp_time(3), temp_time(1)-2000, temp_time(4), temp_time(5));
    params.save_fname = ['movie_save_' tag1];
end

num_iterations = params.num_iterations;
smooth_std = params.smooth_std;
reg_lambda = params.reg_lambda;
high_val_cut_thresh = params.high_val_cut_thresh;
high_val_cut_per_frame = params.high_val_cut_per_frame;
save_all_steps = params.save_all_steps;
save_fname = params.save_fname;
save_dir = params.save_dir;
plot_stuff = params.plot_stuff;

%%

T = size(Y,3);

num_sm_std = size(smooth_std,1);
num_reg_lambda = size(reg_lambda,1);

dsall = cell(num_iterations,1);
Y_reg = Y;

if save_all_steps
    temp_fname = sprintf('%s_pre_moco.h5',save_fname);
    f_save_mov_YS(Y, [save_dir '\' temp_fname], '/mov');
end

%%
if high_val_cut_thresh > 0
    Y_reg = f_mc_cut_high_vals(Y_reg, high_val_cut_thresh, high_val_cut_per_frame);

    if save_all_steps
        temp_fname = sprintf('%s_high_val_cut.h5',save_fname);
        f_save_mov_YS(Y_reg(:,:,1:min(20000,T)), [save_dir '\' temp_fname], '/mov');
    end
end  

for n_iter = 1:num_iterations
    fprintf('Rigid registering iter %d; ', n_iter)

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
        f_save_mov_YS(Y_sm, [save_dir '\' temp_fname], '/mov');
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
        f_save_mov_YS(Y_reg, [save_dir '\' temp_fname_reg], '/mov');
    end
end
mc_out.dsall = dsall;
mc_out.out_frame = mean(Y_reg,3);

mc_out.params = params;
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