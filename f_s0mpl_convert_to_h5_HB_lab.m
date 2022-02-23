function f_s0mpl_convert_to_h5_HB_lab(params)

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
%%
load_type = 2; 
% 1 = Prairie tiffs
% 2 = tiff stack (needs file name)
% 3 = h5 stack (needs file name)

% multiplane data?
params.use_prairie_mpl_tags = 0;
params.mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'}; % 

save_all_steps = 0;
save_indiv_h5info = 0;

% this also saves a trimmed version of movie
params.trim_output_num_frames = 0; %  0 or number of frames to save
params.align_pulse_crop_method = 1;

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

if params.align_pulse_crop_method
    check3 = strfind(params.fname, 'rest');
    if isempty(check3)
        params.align_pulse_crop_method = 1;
    else
        params.align_pulse_crop_method = 2;
    end
end

%%
params_bidi.smooth_std = [1 2 10];%[1 2 2];
params_bidi.fix_range = -50:10;
params_bidi.num_iterations = 1;
params_bidi.plot_stuff = 1;
params_bidi.use_planes = [1 3];

params_moco.image_target = [];
params_moco.num_iterations = 2;
params_moco.plot_stuff = 1;
params_moco.smooth_std = [0.5 0.5 3];

params.params_bidi = params_bidi;
params.params_moco = params_moco;

%%
%load_dir = ['J:\mouse\backup\2018\' file_date '_dLGN\' file_type '-00' file_num];
load_dir = [data_dir '\' params.dset_name];
%load_dir = ['L:\data\Auditory\2018\' file_date '_im\' file_type '-00' file_num];
%load_dir = ['E:\data\V1\' file_date '\' file_type '-00' file_num];

save_dir_movie = [save_dir '\movies'];
save_dir_cuts = [save_dir '\preprocessing'];

save_file_name = params.fname;

disp(save_file_name);

proc_steps = '_cut';

colors1 = parula(num_planes);

params.save_path = save_dir_movie;

if ~isempty(params_moco.im_target_fname)
    moco_init_load = load([save_dir_cuts '\' params_moco.im_target_fname '_h5cutsdata.mat']);
end
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
mat_name = [save_dir_movie '\' save_file_name '_h5cutsdata.mat'];
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
        params.load_path = [load_dir];
        Y_full = bigread3(params.load_path, 1);
    elseif load_type == 3
        params.load_path = [load_dir];
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
    Y_pre_bidi = Y;
    % Y = Y_pre_bidi;
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
        Y{n_pl} = f_bidi_apply_shift(Y{n_pl}, bidi_shifts_all);
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
    Y_pre_moco = Y;
    
    if ~isfield(cuts_data{1}, 'dsall')
        for n_pl = 1:num_planes
            
            if ~isempty(params_moco.im_target_fname)
                if ~isempty(moco_init_load.cuts_data{n_pl}.image_target)
                    params_moco.image_target = moco_init_load.cuts_data{n_pl}.image_target;
                end
            end
            fprintf('%s %s\n', save_file_name, cuts_data{n_pl}.title_tag);
            [~, cuts_data{n_pl}.dsall, cuts_data{n_pl}.image_target] = f_mpl_register2(Y{n_pl}, params_moco);
        end
        save(mat_name, 'params', 'cuts_data');
    end
    
    dsall1 = cell(num_planes, 1);
    for n_pl = 1:num_planes
        dsall1{n_pl} = sum(cat(3,cuts_data{n_pl}.dsall{:}),3);
    end
    dsall1_all = median(cat(3,dsall1{:}),3);
    
    figure;
    subplot(2,1,1); hold on;
    for n_pl = 1:num_planes
        plot(dsall1{n_pl}(:,1), 'color', colors1(n_pl,:));
    end
    plot(dsall1_all(:,1), 'k');
    subplot(2,1,2); hold on;
    for n_pl = 1:num_planes
        plot(dsall1{n_pl}(:,2), 'color', colors1(n_pl,:));
    end
    plot(dsall1_all(:,2), 'k');
    
    for n_pl = 1:num_planes
        Y{n_pl} = uint16(f_suite2p_reg_apply(Y{n_pl}, dsall1_all));
    end
    
    proc_steps = [proc_steps '_moco'];
    if save_all_steps
        for n_pl = 1:num_planes
            f_save_mov_YS(Y{n_pl}, [save_dir_movie '\' save_file_name cuts_data{n_pl}.title_tag proc_steps '.h5'], '/mov')
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