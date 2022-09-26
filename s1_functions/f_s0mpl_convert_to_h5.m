function f_s0mpl_convert_to_h5(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Workflow
%       Load movie, 3 options
%           1: Prairie tiffs
%           2: Tiff stack
%           3: H5 stack
%       Crop pulses
%       Save as H5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% laod params
num_planes = params.num_planes; % number of planes or 0
do_moco = params.do_moco;
do_bidi = params.do_bidi;

data_dir = params.data_dir;
save_dir = params.save_dir;

params_moco.im_target_fname = [params.im_target_fname];
params_moco.do_nonrigid = params.do_nonrigid;
%%
load_type = 1; 
% 1 = Prairie tiffs
% 2 = tiff stack (needs file name)
% 3 = h5 stack (needs file name)

% multiplane data?
if params.num_planes > 1
    params.use_prairie_mpl_tags = 1;
    params.mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'}; % 
else
    params.use_prairie_mpl_tags = 0;
end
save_all_steps = 0;
save_indiv_h5info = 0;

% this also saves a trimmed version of movie
params.trim_output_num_frames = 0; %  0 or number of frames to save

params.overwrite_moco_rigid = 1;
params.overwrite_moco_nonrigid = 1;

colors1 = parula(params.num_planes);

% type 1
%file_type = 'vmmn';close a
%file_type = 'AAF_asynch';
%file_type = 'A1_freq_grating';
%file_type = 'ammn_2_dplanes';
%save_prefix = 'M105_im3_';
%fname = 'AC_ammn_stim'; %
%file_num = '3';
%file_date = '1_21_22a';
%file_date = '10_2_18';
% % type 2 and 3
% load_file_name = 'rest1_5_9_19.hdf5'; % only for 2 and 3

%data_dir = 'D:\data\AC\M105_1_21_22a_dream';
%data_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_24_21_pt3\';

%params_moco.im_target_fname = 'M105_im1_AC_ammn1_1_21_22a_h5cutsdata.mat';%'A1_cont_0.5_12_4_21a_h5cutsdata.mat';

%save_dir = 'C:\Users\rylab_dataPC\Desktop\Yuriy\DD_data\proc_data';
%save_dir = 'E:\data\Auditory\caiman_out_multiplane';
%save_dir = 'J:\mouse\backup\2018\caiman_out_dLGN';
%save_dir = 'L:\data\Auditory\caiman_out';
%save_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data_dream';
%save_dir = 'C:\Users\ys2605\Desktop\stuff\random_save_path';

%%
% 
% if params.align_pulse_crop_method
%     check3 = strfind(params.fname, 'rest');
%     if isempty(check3)
%         params.align_pulse_crop_method = 1;
%     else
%         params.align_pulse_crop_method = 2;
%     end
% end


%%
params_bidi.smooth_std = [1 2 10];%[1 2 2];
params_bidi.fix_range = -50:10;
params_bidi.num_iterations = 1;
params_bidi.plot_stuff = 0;
params_bidi.use_planes = [1 3];

params_moco.image_target = [];
params_moco.plot_stuff = 0;
params_moco.reg_lambda_base = [1 .2];
params_moco.high_val_cut_thresh = 0.01;

if params.moco_rigid_method == 1 % regular multiplane
    params_moco.num_iterations = 2;
    
    params_moco.smooth_std = [0.5 0.5 6;...
                              0.5 0.5 3;... % was 3 for missmatch
                              0.5 0.5 1;...
                              0.5 0.5 0.5];

    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              2 .5];
                          
elseif params.moco_rigid_method == 2 % regular missmatch 30 hz
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0 0 0.5];
                          
    params_moco.reg_lambda = [0 .2;... % 1
                              2 .2;...
                              2 .5;...
                              2 .5];
elseif params.moco_rigid_method == 21 % regular missmatch 30 hz
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0 0 0.5];
                          
    params_moco.reg_lambda = [1 .2;... % 1
                              2 .2;...
                              2 .5;...
                              2 .5];
                          
elseif params.moco_rigid_method == 22 % noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 2;...
                              0.5 0.5 2];
                          
    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              2 .5];
                          
elseif params.moco_rigid_method == 23 % even more noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 5; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 .5;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              4 1;...
                              4 1];
elseif params.moco_rigid_method == 24 % even more noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 2; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 .5;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [5 .5;...
                              2 .2;...
                              2 .5;...
                              4 1;...
                              4 1];
                          
elseif params.moco_rigid_method == 25 % noisy missmatch 30 hz, 
    
    params_moco.num_iterations = 6; % 4 was for mmn data works with 30hz noisy data

    params_moco.smooth_std = [0.5 0.5 6;...
                              0.5 0.5 3;... % was 3 for missmatch
                              0.5 0.5 2;...
                              0.5 0.5 1;...
                              0.5 0.5 0;...
                              0.5 0.5 0];
                          
    params_moco.reg_lambda = [.1 .01;...
                              .1 .01;...
                              .1 .01;...
                              .1 .01];
                                                
elseif params.moco_rigid_method == 3 % multiplane super noisy; dream/chrmine
    
    params_moco.num_iterations = 4; % 

    params_moco.smooth_std = [0.5 0.5 12;...
                              0.5 0.5 6;... % was 3 for missmatch
                              0.5 0.5 3;...
                              0.5 0.5 3];
                          
    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              2 .5];
                          
elseif params.moco_rigid_method == 32 % even more super noisy; dream/chrmine
    
    params_moco.num_iterations = 3; % 

    params_moco.smooth_std = [1 1 12;...
                              1 1 7;... 
                              1 1 5;...
                              0.5 0.5 3];
                          
    params_moco.reg_lambda = [1 .2;...
                              2 .2;...
                              2 .5;...
                              2 .5];
end
params_moco.medfilt = 0;
                      
% % for better snr
% params_moco.num_iterations = 5;
% params_moco.smooth_std = [0.5 0.5 6;...        
%                           0.5 0.5 3;...
%                           0 0 2;];       % each row corresponds to smooth iterations

% params_moco.num_iterations = 8;
% params_moco.smooth_std = [0.5 0.5 12;...
%                           0.5 0.5 12;...
%                           0.5 0.5 9;...
%                           0.5 0.5 9;...
%                           0.5 0.5 7;...
%                           0.5 0.5 7;...
%                           0.5 0.5 5;...
%                           0.5 0.5 5];       % each row corresponds to smooth iterations
% params_moco.reg_lambda = [1e0; 1e0; 5e0; 5e0; 5e0; 5e0; 1e1; 1e1];
%                       

if params.moco_nonrigid_method == 1
    params_moco.nonrigid_smooth_std = [0.5 0.5 6];
    params_moco.nonrigid_reg_lambda = [2 .5];
    params_moco.nonrigid_block_size = 60;
    params_moco.nonrigid_block_overlap = 40;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3];
    
elseif params.moco_nonrigid_method == 2
    params_moco.nonrigid_smooth_std = [0.5 0.5 1];
    params_moco.nonrigid_reg_lambda = [2 .5];
    params_moco.nonrigid_block_size = 50;
    params_moco.nonrigid_block_overlap = 30;
    params_moco.nonrigid_block_smooth = [0.5 0.5 1];

elseif params.moco_nonrigid_method == 3
    params_moco.nonrigid_smooth_std = [0.5 0.5 3];
    params_moco.nonrigid_reg_lambda = [2 .5];
    params_moco.nonrigid_block_size = 40;
    params_moco.nonrigid_block_overlap = 30;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3]; % [0 0 025]
elseif params.moco_nonrigid_method == 4
    params_moco.nonrigid_smooth_std = [0.5 0.5 3];
    params_moco.nonrigid_reg_lambda = [2 .5];
    params_moco.nonrigid_block_size = 30;
    params_moco.nonrigid_block_overlap = 15;
    params_moco.nonrigid_block_smooth = [0.5 0.5 3]; % [0 0 025]
    
end


%%
%load_dir = ['J:\mouse\backup\2018\' file_date '_dLGN\' file_type '-00' file_num];
load_dir = [data_dir '\' params.dset_name(1:end-1) '-00' params.dset_name(end) '\'];
%load_dir = ['L:\data\Auditory\2018\' file_date '_im\' file_type '-00' file_num];
%load_dir = ['E:\data\V1\' file_date '\' file_type '-00' file_num];

save_dir_movie = [save_dir '\movies'];
save_dir_cuts = [save_dir '\preprocessing'];

save_file_name = params.fname;

disp(save_file_name);
fprintf('Pulse corp method = %d\n', params.align_pulse_crop_method);
fprintf('Moco rigid method = %d\n', params.moco_rigid_method);
if params.do_nonrigid
    fprintf('Moco nonrigind method = %d\n', params.moco_nonrigid_method);
end

proc_steps = '_cut';

params.save_path = save_dir_movie;

if ~isempty(params_moco.im_target_fname)
    moco_init_load = load([save_dir_cuts '\' params_moco.im_target_fname '_h5cutsdata.mat']);
end

params_moco.save_fname = params.fname;
params_moco.save_dir = save_dir_movie;
params.params_bidi = params_bidi;
params.params_moco = params_moco;

%%
if ~exist(save_dir_movie, 'dir')
    mkdir(save_dir_movie)
end
if ~exist(save_dir_cuts, 'dir')
    mkdir(save_dir_cuts)
end
if ~exist([save_dir_movie '\ave_proj'], 'dir')
    mkdir([save_dir_movie '\ave_proj'])
end

addpath([pwd '\general_functions']);

%% load cuts data
mat_name = [save_dir_cuts '\' save_file_name '_h5cutsdata.mat'];
if exist(mat_name, 'file')
    load_data = load(mat_name);
    cuts_data = load_data.cuts_data;
else
    cuts_data = cell(num_planes,1);
end

%% load
Y = cell(num_planes,1);

if params.use_prairie_mpl_tags
    params.load_path = load_dir;
    for n_pl = 1:num_planes
        Y{n_pl} = f_collect_prairie_tiffs4(params.load_path, params.mpl_tags{n_pl});
    end
else
    if load_type == 1
        params.load_path = load_dir;
        Y_full = f_collect_prairie_tiffs4(params.load_path, 'Ch2');
    elseif load_type == 2
        params.load_path = [load_dir, '\',  load_file_name];
        Y_full = bigread3(params.load_path, 1);
    elseif load_type == 3
        params.load_path = [load_dir, '\',  load_file_name];
        Y_full = h5read(params.load_path, '/mov');
    end
    if num_planes > 1
        last_time = size(Y_full,3);
        params.ave_trace_full = squeeze(mean(mean(Y_full, 1),2));
        figure; plot(params.ave_trace_full)
        title('Full ave trace');
        for n_pl = 1:num_planes
            ind_mpl = n_pl:num_planes:last_time;
            Y{n_pl} = Y_full(:,:,ind_mpl);
        end
    else
        Y{1} = Y_full;
    end
    clear Y_full;
end

%% compute cuts
if ~isfield(cuts_data{1}, 'vid_cuts_trace')
    cuts_data = cell(num_planes,1);
    [d1, d2, T] = size(Y{1});
    vid_cuts_trace_all = true(T,1);
    for n_pl = 1:num_planes
        cuts_data{n_pl} = params;
        if num_planes>1
            cuts_data{n_pl}.title_tag = sprintf('_mpl%d_pl%d', num_planes, n_pl);
        else
            cuts_data{n_pl}.title_tag = '';
        end
        cuts_data{n_pl}.ave_trace = squeeze(mean(mean(Y{n_pl}, 1),2));
        cuts_data{n_pl} = f_s0_compute_align_cuts(cuts_data{n_pl});
        vid_cuts_trace_all = vid_cuts_trace_all.*cuts_data{n_pl}.vid_cuts_trace;
    end
    for n_pl = 1:num_planes
        cuts_data{n_pl}.vid_cuts_trace = vid_cuts_trace_all;
    end
    save(mat_name, 'params', 'cuts_data');
end

%% apply cuts
for n_pl = 1:num_planes
    Y{n_pl}(:,:,~cuts_data{n_pl}.vid_cuts_trace) = [];
    if save_all_steps
        f_save_mov_YS(Y{n_pl}, [save_dir_movie '\' save_file_name cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov');
    end
end

%% bidi fix
if do_bidi
    Y_pre_corr = Y;
    % Y = Y_pre_corr;
    if ~isfield(cuts_data{1}, 'bidi_out')
        % compute
        for n_pl = 1:num_planes
            fprintf('%s %s\n', save_file_name, cuts_data{n_pl}.title_tag);
            params_bidi.title_tag = cuts_data{n_pl}.title_tag;
            [~, cuts_data{n_pl}.bidi_out] = f_fix_bidi_shifts3(Y{n_pl}, params_bidi);
        end
        save(mat_name, 'params', 'cuts_data');
    end
    
    bidi_shifts_all = cell(num_planes,1);
    for n_pl = 1:num_planes
        bidi_shifts_all{n_pl} = sum(cuts_data{n_pl}.bidi_out.best_shifts,2);
    end
    bidi_shifts_all = cat(2,bidi_shifts_all{:});
    % only use top 3
    if isfield(params_bidi, 'use_planes')
        mean_tag = num2str(params_bidi.use_planes(1):min([params_bidi.use_planes(2) num_planes]));
        mean_bidi_shifts = round(mean(bidi_shifts_all(:,params_bidi.use_planes(1):min([params_bidi.use_planes(2) num_planes])),2));
    else
        mean_tag = 'all';
        mean_bidi_shifts = round(mean(bidi_shifts_all(:,1:min([3 num_planes])),2));
    end
    
    figure; hold on;
    for n_pl = 1:num_planes
        plot(bidi_shifts_all(:,n_pl), 'color', colors1(n_pl, :))
    end
    plot(mean_bidi_shifts, 'k');
    title(['computed bidi shifts all planes; black-average pl ' mean_tag])
    

    % apply
    for n_pl = 1:num_planes
        [d1, d2, T] = size(Y{n_pl});
        Y_temp = single(Y{n_pl});
        if cuts_data{n_pl}.bidi_out.params.do_interp
            deg_per_fov = 180 * cuts_data{n_pl}.bidi_out.params.laser_open_frac;
            y0 = 1:d1;
            z0 = 1:T;
            deg0 = linspace(-deg_per_fov/2, deg_per_fov/2, d2);
            x0 = sin(deg0/360*2*pi);
            x0n = x0 - min(x0);
            x0n = x0n/max(x0n)*(d2-1)+1;
            x_coords = 1:d2;
            [X_corr, Y0, Z0] = meshgrid(x_coords, y0, z0);
            [X_real, ~, ~] = meshgrid(x0n, y0, z0);
            Y_temp = interp3(X_corr, Y0, Z0, Y_temp, X_real, Y0, Z0);
        end
        Y_temp = f_bidi_apply_shift(Y_temp, bidi_shifts_all);
        if cuts_data{n_pl}.bidi_out.params.do_interp
            Y_temp = interp3(X_real, Y0, Z0, Y_temp, X_corr, Y0, Z0);
        end
        Y{n_pl} = uint16(Y_temp);
    end
    proc_steps = [proc_steps '_bidi'];
    
    if save_all_steps
        for n_pl = 1:num_planes
            f_save_mov_YS(Y{n_pl}, [save_dir_movie '\' save_file_name cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov');
        end
    end
end

%% moco
if do_moco
    Y_pre_corr = Y;
    %Y = Y_pre_corr;
    
    % correct movie to itself
    for n_pl = 1:num_planes
        if ~isfield(cuts_data{n_pl}, 'dsall') || params.overwrite_moco_rigid
            fprintf('%s %s\n', save_file_name, cuts_data{n_pl}.title_tag);
            [~, mc_out] = f_mc_rigid(Y{n_pl}, params_moco);
            cuts_data{n_pl}.dsall = mc_out.dsall;
        end
    end
    
    save(mat_name, 'params', 'cuts_data');
    
    % add all coorection and use median actoss planes
    [~, dsall1_all, dsall1_all_mf] = f_mc_dsall_proc(cuts_data);
    
    if params_moco.medfilt
        dsall1_use = dsall1_all_mf;
    else
        dsall1_use = dsall1_all;
    end
    f_mc_plot_cuts_data(cuts_data, save_file_name);
    
    %Y2 = Y_pre_moco;
    for n_pl = 1:num_planes
        Y{n_pl} = uint16(f_suite2p_reg_apply(Y{n_pl}, dsall1_use));
        %Y2{n_pl} = uint16(f_suite2p_reg_apply(Y_pre_moco{n_pl}, dsall1_all_r));
    end
    
    % apply global offset to while movie
    for n_pl = 1:num_planes
        Y_temp = single(Y{n_pl});
        cuts_data{n_pl}.image_target = mean(Y_temp,3);
        cuts_data{n_pl}.image_target_std = std(Y_temp,0,3);
    end

    ds_base_all = zeros(num_planes, 2);
    for n_pl = 1:num_planes
        cuts_data{n_pl}.ds_base = [0 0];
        if ~isempty(params_moco.im_target_fname)
            if ~isempty(moco_init_load.cuts_data{n_pl}.image_target)
                cuts_data{n_pl}.image_target_external = moco_init_load.cuts_data{n_pl}.image_target;
                cuts_data{n_pl}.ds_base = f_suite2p_reg_compute(cuts_data{n_pl}.image_target, cuts_data{n_pl}.image_target_external, params_moco.reg_lambda_base);
            end
        end
        ds_base_all(n_pl, :) = cuts_data{n_pl}.ds_base;
    end
            
    save(mat_name, 'params', 'cuts_data');

    if ~isempty(params_moco.im_target_fname)
        figure; plot(ds_base_all); title('correction to external input database');
    end
    
    dsall1_use2 = ones(size(dsall1_use));
    for n_pl = 1:num_planes
        Y{n_pl} = uint16(f_suite2p_reg_apply(Y{n_pl}, dsall1_use2.*ds_base_all(n_pl, :)));
        %Y2{n_pl} = uint16(f_suite2p_reg_apply(Y_pre_moco{n_pl}, dsall1_all_r));
    end
     
    if params.moco_zero_edge
        for n_pl = 1:num_planes
            Y{n_pl} = f_mc_zero_edges(Y{n_pl}, dsall1_use, ds_base_all(n_pl, :));
        end
    end
    
    proc_steps = [proc_steps '_moco'];
    if save_all_steps
        for n_pl = 1:num_planes
            f_save_mov_YS(Y{n_pl}(:,:,1:min(25000, size(Y{n_pl},3))), [save_dir_movie '\' save_file_name cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov')
        end
    end
    
    if params.do_nonrigid
        Y_pre_corr = Y;
        %Y = Y_pre_corr;
        
        for n_pl = 1:num_planes
            if ~isfield(cuts_data{n_pl}, 'nr_corr_data') || params.overwrite_moco_nonrigid
                [~, mc_out2] = f_mc_nonrigid(Y{n_pl}, params_moco);
                cuts_data{n_pl}.nr_corr_data = mc_out2.nr_corr;
            end
        end

        f_mc_plot_nrcuts_data(cuts_data, save_file_name)
        
        for n_pl = 1:num_planes
            Y{n_pl} = f_mc_apply_nonrigid_corr(Y{n_pl}, cuts_data{n_pl}.nr_corr_data);
        end
        proc_steps = [proc_steps '_nonrigid'];
        if 1%save_all_steps
            for n_pl = 1:num_planes
                f_save_mov_YS(Y{n_pl}(:,:,1:min(25000, size(Y{n_pl},3))), [save_dir_movie '\' save_file_name cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov')
            end
        end
    end
    
end

%% save
params.params_bidi = params_bidi;
params.params_moco = params_moco;

for n_pl = 1:num_planes
    params.save_mov_path = [save_dir_movie '\' save_file_name cuts_data{n_pl}.title_tag '.h5'];
    params.cuts_data = cuts_data{n_pl};
    
    f_save_mov_YS(Y{n_pl}, params.save_mov_path, '/mov');
    
    if params.trim_output_num_frames
        params.save_mov_path_trim = [save_dir_movie '\' save_file_name cuts_data{n_pl}.title_tag proc_steps '_trim.hdf5'];
        f_save_mov_YS(Y{n_pl}(:,:,1:round(params.trim_output_num_frames)), params.save_mov_path_trim, '/mov');  
    end
    
    if save_indiv_h5info
        save([save_dir '\' save_file_name cuts_data{n_pl}.title_tag '_h5cutsinfo.mat'], 'params');
    end
    
    tmp_fig = figure; imagesc(squeeze(mean(Y{n_pl},3)));
    title([save_file_name ' Ave prjection ' cuts_data{n_pl}.title_tag], 'Interpreter', 'none');
    axis tight equal;
    saveas(tmp_fig,[save_dir_movie '\ave_proj\' save_file_name cuts_data{n_pl}.title_tag '_ave_proj.fig']);
end

disp('Done')
%% analysis

end