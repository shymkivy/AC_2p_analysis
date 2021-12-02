function [Y_reg, dsall] = f_mpl_register(Y, params)
% all inputs and outputs should be 1d cells

if ~exist('params', 'var')
    params = struct();
end

if isfield(params, 'image_target')
    image_target = params.image_target;
else
    image_target = [];
end

if isfield(params, 'save_smooth')
    save_smooth = params.save_smooth;
else
    save_smooth = 0;
end

if isfield(params, 'save_fname')
    save_fname = params.save_fname;
else
    temp_time = clock;
    tag1 = sprintf('%d_%d_%d_%dh_%dm',temp_time(2), temp_time(3), temp_time(1)-2000, temp_time(4), temp_time(5));
    save_fname = ['movie_save_' tag1];
end

%%
num_planes = numel(Y);
dsall = cell(num_planes,1);

for n_pl = 1:num_planes       
    % smooth movie
    smooth_std = [0.5 0.5 3];
    %smooth_std = [0 0 3];
    Y_sm = f_smooth_movie(Y{n_pl}, smooth_std);

    %%
    if save_smooth
        temp_fname_sm = sprintf('%s_pl%d_sm.h5',save_fname, n_pl);
        f_save_mov_YS(Y_sm, [save_dir_movie '\' temp_fname_sm], '/mov');
    end
    %%
    tic;
    [~, dsall{n_pl}] = f_suite2p_register_YS(Y_sm, image_target{n_pl});
    toc
    clear Y_sm;

%         figure; plot(dsall{n_pl})
%         
%         temp_fname_reg = sprintf('%s_%s_pl%d_reg.h5',fname, file_date, n_pl);
%         f_save_mov_YS(Y_reg, [save_dir_movie '\' temp_fname_reg], '/mov');

end

%% apply corr to raw files
Y_reg = cell(multiplane, 1);
for n_pl = 1:multiplane
    tic;
    Y_reg{n_pl} = f_suite2p_apply_reg_YS(Y{n_pl}, dsall{n_pl});
    toc
end


end