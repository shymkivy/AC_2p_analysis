function [Y_reg, mc_out] = f_mc_nonrigid(Y, params)

if ~isfield(params, 'nonrigid_block_size'); params.nonrigid_block_size = 60; end
if ~isfield(params, 'nonrigid_block_overlap'); params.nonrigid_block_overlap = 40; end
if ~isfield(params, 'nonrigid_block_smooth'); params.nonrigid_block_smooth = [1, 1, 1]; end
if ~isfield(params, 'nonrigid_smooth_std'); params.nonrigid_smooth_std = [.5 .5 3]; end
if ~isfield(params, 'nonrigid_reg_lambda'); params.nonrigid_reg_lambda = [2 .5]; end
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

nr_smooth_std = params.nonrigid_smooth_std;
high_val_cut_thresh = params.high_val_cut_thresh;
high_val_cut_per_frame = params.high_val_cut_per_frame;
save_all_steps = params.save_all_steps;
save_fname = params.save_fname;
save_dir = params.save_dir;
plot_stuff = params.plot_stuff;

%%
T = size(Y,3);
mc_out.out_frame = mean(Y,3);

params_nrcorr.target_frame = mc_out.out_frame;
params_nrcorr.reg_lambda = params.nonrigid_reg_lambda;
params_nrcorr.nonrigid_block_size = params.nonrigid_block_size;
params_nrcorr.nonrigid_block_overlap = params.nonrigid_block_overlap;
params_nrcorr.nonrigid_block_smooth = params.nonrigid_block_smooth;

%%
fprintf('Nonrigid registering; ')

Y_reg = Y;

if high_val_cut_thresh > 0
    tic;
    Y_reg = f_mc_cut_high_vals(Y_reg, high_val_cut_thresh, high_val_cut_per_frame);
    fprintf('high val cut %.1f%% duration=%.1fsec; ', high_val_cut_thresh*100, toc);
    if save_all_steps
        temp_fname = sprintf('%s_nrcorr_high_val_cut.h5',save_fname);
        f_save_mov_YS(Y_reg(:,:,1:min(20000,T)), [save_dir '\' temp_fname], '/mov');
    end
end  

if sum(nr_smooth_std>0)
    tic;
    Y_reg = f_smooth_movie(Y_reg, nr_smooth_std); % uses ram
    %Y_sm = f_smooth_movie2(Y, smooth_std); % uses GPU
    fprintf('smooth [%.1f %.1f %.1f] duration=%.1fsec; ', nr_smooth_std(1), nr_smooth_std(2), nr_smooth_std(3), toc);
    
    if save_all_steps
        temp_fname = sprintf('%s_nrcorr_sm.h5',save_fname);
        f_save_mov_YS(Y_reg(:,:,1:min(20000,T)), [save_dir '\' temp_fname], '/mov');
    end
end
fprintf('\n')

nr_corr = f_mc_compute_nonrigid_motion(Y_reg, params_nrcorr);

if save_all_steps
    Y_reg = f_mc_apply_nonrigid_corr(Y, nr_corr);
    temp_fname_reg = sprintf('%s_nrcorr_post.h5',save_fname);
    f_save_mov_YS(Y_reg(:,:,1:min(20000,T)), [save_dir '\' temp_fname_reg], '/mov');
end

mc_out.nr_corr = nr_corr;


end