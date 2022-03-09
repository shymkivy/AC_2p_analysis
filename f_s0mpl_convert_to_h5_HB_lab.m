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
params_moco.im_target_fname = [params.im_target_fname];
%%

% multiplane data?
params.mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'}; % 

save_all_steps = 0;
save_indiv_h5info = 0;

% this also saves a trimmed version of movie
params.trim_output_num_frames = 0; %  0 or number of frames to save
params.align_pulse_crop_method = 0; % 0 = no cuts, 1 = auto cut; 2 = manual cuts

%%

if params.align_pulse_crop_method
    check3 = strfind(params.fname, 'rest');
    if isempty(check3)
        params.align_pulse_crop_method = 1;
    else
        params.align_pulse_crop_method = 0;
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
load_dir = params.data_dir;
%load_dir = ['L:\data\Auditory\2018\' file_date '_im\' file_type '-00' file_num];
%load_dir = ['E:\data\V1\' file_date '\' file_type '-00' file_num];

[~, f_name_core, f_ext] = fileparts(params.fname);

save_dir_movie = params.save_dir;
save_dir_cuts = params.save_dir;

save_file_name = f_name_core;

disp(save_file_name);

proc_steps = '_cut';

colors1 = parula(num_planes);

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
    if params.load_type == 1
        params.load_path = [load_dir '\' params.fname];
        Y_full = f_collect_prairie_tiffs4(params.load_path, 'Ch2');
    else
        params.load_path = [load_dir '\' params.fname];
        if or(strcmpi(f_ext, '.h5'), strcmpi(f_ext, '.hdf5'))
            Y_full = h5read(params.load_path, '/mov');
        elseif or(strcmpi(f_ext, '.tif'), strcmpi(f_ext, '.tiff'))
            Y_full = bigread3(params.load_path, 1);
        end
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
%% try pmt nouse fix?

[d1, d2 , T] = size(Y{1});

figure; imagesc(mean(Y{1},3)')

Y_flat = reshape(Y{1}, [d1*d2*T, 1]);

Y_flat2 = (Y_flat - mean(Y_flat))/std(Y_flat);

figure; plot(Y_flat)

fs = d1*d1;

window = 1024;
noverlap = 120;

[s,f,t] = spectrogram(Y_flat2,window,noverlap,[], fs);

figure;
imagesc(t, f, real(s))


%%

shifts_range = 200;

frame1 = Y{1}(:,:,1);
frame1n = (frame1 - mean(frame1(:)))/std(frame1(:));
frame1n2 = frame1n;
frame1n2(frame1n>0.5) = 0.5;
frame1n2(frame1n<-0.5) = -0.5;
for n_fr = 2:T
    frame_temp = Y{1}(:,:,n_fr);
    frame_tempn = (frame_temp - mean(frame_temp(:)))/std(frame_temp(:));
    frame_tempn2 = frame_tempn;
    frame_tempn2(frame_tempn>0.5) = 0.5;
    frame_tempn2(frame_tempn<-0.5) = -0.5;

    corr_vals = zeros(shifts_range+1,1);
    
    corr_vals(1) = mean(mean((frame1n2 .* frame_tempn2)));
    for sh1 = 1:shifts_range
    	corr_vals(sh1+1) = mean(mean((frame1n2 .* circshift(frame_tempn2, [0 sh1]))));
    end
    
end

figure; plot(corr_vals)

figure; imagesc(frame1n2)

fft_all = zeros(d1, d2, T);
for n_fr = 1:T
    frame_temp = Y{1}(:,:,n_fr);
    frame_temp = (frame_temp - mean(frame_temp(:)))/std(frame_temp(:));
    fft_all(:,:,n_fr) = (fft2(frame_temp));
    
end

figure; imagesc(ifft2(mean(fft_all,3)))

fft2()

frame1 = Y{1}(:,:,4844);

frame2 = Y{1}(:,:,4855);

frame11 = frame1;
frame11(frame1<3140) = 3140;
frame11(frame1>3160) = 3160;
frame11 = (frame11 - mean(frame11(:)))/std(frame11(:));

frame22 = frame2;
frame22(frame2<3140) = 3140;
frame22(frame2>3160) = 3160;
frame22 = (frame22 - mean(frame22(:)))/std(frame22(:));

frame111 = frame1;
frame111 = (frame111 - mean(frame111(:)))/std(frame111(:));
frame222 = frame2;
frame222 = (frame222 - mean(frame222(:)))/std(frame222(:));

num_shifts = 300;
per1 = 100;
out1 = zeros(num_shifts,1);
for n_px = 1:num_shifts
    out1(n_px) = mean(mean(frame111(100:(100+per1),:) .* frame222((n_px+50):(n_px+50+per1),:)));
end

figure; plot(out1)

frame1_ff1 = fft2(frame111);

frame1_ff2 = fft2(frame222);

frame_iff = ifft2(frame1_ff1.*frame1_ff2);

figure; imagesc(frame_iff)

figure; plot(sum(frame_iff))

figure; imagesc(abs(fftshift(frame1_ff)))


figure; imagesc(ifft2(fft2(frame1)))

ifft2(frame1_ff)

figure; imagesc(frame1')
caxis([3100 3200])

figure; imagesc(frame2')
caxis([3100 3200])

figure; plot(frame2(:))

frames = Y{1}(:,:,2255:2260);

frames2 = (frames - mean(frames(:)))/std(frames(:));

frame2 = (frame1 - mean(frame1(:)))/std(frame1(:));

fs = 512;
window = 128;
noverlap = 120;

figure;
spectrogram(frames2(:),window,noverlap,[], fs, 'yaxis');

fftout = fftshift(fft(frame1(:)));

figure; plot(abs(fftout))


[num,den] = iirnotch(14,1);


[s,f,t] = spectrogram(Y_flat2,window,noverlap,[], fs);


pxx = pwelch(frame2(:),window, noverlap);

figure; plot(pxx)

%% bidi fix
if params.do_bidi
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
if params.do_moco
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