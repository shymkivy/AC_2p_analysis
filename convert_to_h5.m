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
%%
clear;
%close all;

%%
load_type = 1; 
% 1 = Prairie tiffs
% 2 = tiff stack (needs file name)
% 3 = h5 stack (needs file name)

% multiplane data?
multiplane = 5; % number of planes or 0

params.auto_align_pulse_crop = 1;
% this also saves a trimmed version of movie
params.trim_output_num_frames = 0; %  0 or number of frames to save

data_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\1_30_21_im_stim\';

% type 1
%file_type = 'vmmn';
%file_type = 'AAF_asynch';
%file_type = 'A1_freq_grating';
%file_type = 'ammn_2_dplanes';
file_name = 'AC_ammn_5plt_1plx'; % 
file_num = '2';
file_date = '1_30_21';
%file_date = '10_2_18';
%
%load_dir = ['J:\mouse\backup\2018\' file_date '_dLGN\' file_type '-00' file_num];
load_dir = [data_dir '\' file_name '-00' file_num];
%load_dir = ['L:\data\Auditory\2018\' file_date '_im\' file_type '-00' file_num];
%load_dir = ['E:\data\V1\' file_date '\' file_type '-00' file_num];

params.use_prairie_mpl_tags = 1;
params.mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'}; % 

% % type 2 and 3
% load_file_name = 'rest1_5_9_19.hdf5'; % only for 2 and 3

%save_dir = 'C:\Users\rylab_dataPC\Desktop\Yuriy\DD_data\proc_data';
%save_dir = 'E:\data\Auditory\caiman_out_multiplane';
%save_dir = 'J:\mouse\backup\2018\caiman_out_dLGN';
%save_dir = 'L:\data\Auditory\caiman_out';
save_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';
%save_dir = 'C:\Users\ys2605\Desktop\stuff\random_save_path';

save_dir_movie = [save_dir '\movies'];

save_file_name = [file_name file_num '_' file_date];
disp(save_file_name);


%%
if ~exist(save_dir_movie, 'dir')
    mkdir(save_dir_movie)
end
if ~exist([save_dir_movie '\ave_proj'], 'dir')
    mkdir([save_dir_movie '\ave_proj'])
end

addpath([pwd '\general_functions']);

%% load

    
    
if ~params.use_prairie_mpl_tags
    if load_type == 1
        params.load_path = load_dir;
        Y = f_collect_prairie_tiffs4(params.load_path, 'Ch2');
    elseif load_type == 2
        params.load_path = [load_dir, '\',  load_file_name];
        Y = bigread3(params.load_path, 1);
    elseif load_type == 3
        params.load_path = [load_dir, '\',  load_file_name];
        Y = h5read(params.load_path, '/mov');
    end
end

%%
if multiplane
    if params.use_prairie_mpl_tags
        params.load_path = load_dir;
        Y = cell(multiplane,1);
        for n_pl = 1:multiplane
            Y{n_pl} = f_collect_prairie_tiffs4(params.load_path, params.mpl_tags{n_pl});
        end
    else
        last_time = size(Y,3);
        params.ave_trace_full = squeeze(mean(mean(Y, 1),2));
        figure; plot(params.ave_trace_full)
        title('Full ave trace');
    end
    
    Y_full = Y;

    for n_pl = 1:multiplane
        if params.use_prairie_mpl_tags
            Y = Y_full{n_pl};
        else
            ind_mpl = n_pl:multiplane:last_time;
            Y = Y_full(:,:,ind_mpl);
        end
        
        params.ave_trace = squeeze(mean(mean(Y, 1),2));

        params = if_compute_align_cuts(params);
        
        Y(:,:,~logical(params.vid_cuts_trace)) = [];
        
        params.save_mov_path = [save_dir_movie '\' save_file_name '_mpl' num2str(n_pl) '_cut.hdf5'];
        f_save_mov_YS(Y, params.save_mov_path, '/mov');
        if params.trim_output_num_frames
            params.save_mov_path_trim = [save_dir_movie '\' save_file_name '_mpl' num2str(n_pl) '_cut_trim.hdf5'];
            f_save_mov_YS(Y(:,:,1:round(params.trim_output_num_frames)), params.save_mov_path_trim, '/mov');  
        end
        save([save_dir '\' save_file_name '_mpl' num2str(n_pl) '_h5cutsinfo.mat'], 'params');

        tmp_fig = figure; imagesc(squeeze(mean(Y,3))); title([save_file_name ' Ave prjection multiplane pl' num2str(n_pl)], 'Interpreter', 'none'); axis tight equal;
        saveas(tmp_fig,[save_dir_movie '\ave_proj\' save_file_name '_mpl' num2str(n_pl) '_ave_proj']);
    end
else
    %% process
    % mean trace
    params.ave_trace = squeeze(mean(mean(Y, 1),2));

    [params] = if_compute_align_cuts(params);

    % crop movie
    Y(:,:,~logical(params.vid_cuts_trace)) = [];

    %% save h5 file
    params.save_mov_path = [save_dir_movie '\' save_file_name '_cut.h5'];

    f_save_mov_YS(Y, params.save_mov_path, '/mov');
    if params.trim_output_num_frames
        params.save_mov_path_trim = [save_dir_movie '\' save_file_name '_cut_trim.hdf5'];
        f_save_mov_YS(Y(:,:,1:round(params.trim_output_num_frames)), params.save_mov_path_trim, '/mov');
    end
    save([save_dir '\' save_file_name '_h5cutsinfo.mat'], 'params');
    tmp_fig = figure; imagesc(squeeze(mean(Y,3))); title([save_file_name ' Ave prjection'], 'Interpreter', 'none'); axis equal tight;
    saveas(tmp_fig,[save_dir_movie '\' save_file_name '_ave_proj']);
end

disp('Done')
%% analysis

%% functions
function [params] = if_compute_align_cuts(params)
    if ~isfield(params, 'auto_align_pulse_crop')
        params.auto_align_pulse_crop = 1; % default
    end
    ave_trace = params.ave_trace;

    min_trace = min(ave_trace);
    max_trace = max(ave_trace);
    norm_ave_trace = (ave_trace - min_trace)/(max_trace-min_trace);
    
    if params.auto_align_pulse_crop
        thresh = 0.5;
        pulse_buff = 30; % frames
        
        pulse_trace = (norm_ave_trace);
        pulse_trace(pulse_trace<thresh) = 0;
        pulse_trace(pulse_trace>thresh) = 1;
       
        pulse_on = find(diff(pulse_trace)>0)+1;
        pulse_off = find(diff(pulse_trace)<0)+1;
        
        % quality check if light turns off in beginning 
        if and(pulse_trace(1) == 1,numel(pulse_off)>numel(pulse_on))
            pulse_off(1) = [];
        end
        
        num_frag = round(numel(pulse_on)-1);
        
        vid_cuts = zeros(num_frag,2);
        
        n_pulse = 1;
        vid_cuts_trace = zeros(size(norm_ave_trace));
        for n_frag = 1:num_frag
            vid_cuts(n_frag,1) = pulse_off(n_pulse) + pulse_buff;
            n_pulse = n_pulse+1;
            vid_cuts(n_frag,2) = pulse_on(n_pulse) - pulse_buff;
            vid_cuts_trace(vid_cuts(n_frag,1):vid_cuts(n_frag,2)) = 1;
        end

        
        figure; plot(norm_ave_trace);
        hold on; plot(vid_cuts_trace);
        axis tight;
        title(sprintf('Automatic %d frag selected', num_frag));
    else
        f1 = figure;
        plot(norm_ave_trace);
        axis tight;
        title('how many fragments?');
        num_frag = input('how many fragments? (int):');
        vid_cuts = zeros(num_frag,2);
        vid_cuts_trace = zeros(size(norm_ave_trace));
        for n_frag = 1:num_frag
            title(sprintf('Select fragment %d/%d (2 clicks)', n_frag,num_frag));
            [temp_cuts, ~] = ginput(2);
            vid_cuts(n_frag,:) = round(temp_cuts);
            if vid_cuts(n_frag,1) < 1
                vid_cuts(n_frag,1) = 1;
            end
            if vid_cuts(n_frag,2) > numel(ave_trace)
                vid_cuts(n_frag,2) = numel(ave_trace);
            end
            vid_cuts_trace(vid_cuts(n_frag,1):vid_cuts(n_frag,2)) = 1;
            plot(norm_ave_trace);
            hold on;
            plot(vid_cuts_trace);
            hold off;
            axis tight;
        end
        close(f1)
        
        figure;plot(norm_ave_trace);
        hold on; plot(vid_cuts_trace);
        axis tight;
        title(sprintf('Manual %d frag selected', num_frag));
    end
    
    params.vid_cuts = vid_cuts;
    params.vid_cuts_trace = vid_cuts_trace;
end
