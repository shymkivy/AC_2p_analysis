close all
clear

%%
addpath([pwd '\s1_functions']);
addpath([pwd '\general_functions'])

%%
%ops.file_dir = 'F:\AC_data\caiman_data_dream3\preprocessing';
%ops.file_dir = 'F:\AC_data\caiman_data_echo\preprocessing';

params.dset_table_fpath = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\AC_data_list_all.xlsx';
params.data_dir = 'F:\AC_data\caiman_data_dream3';
%params.data_dir = 'D:\data\caiman_data_dream';

% these have to be same as column names in excel
params.limit.dset_name =        '';
params.limit.experiment =       'dream';
params.limit.mouse_id =         'M166';
params.limit.mouse_tag =        '6_20_22_pt2';
params.limit.dset_name =        '';

roi_half_size = 6;

%%
AC_data = f_s0_parse_tab_data(params);

%%
mouse_id_all = unique(AC_data.mouse_id, 'stable');

% work plane by plane

num_dsets = numel(AC_data.mouse_id);
num_planes = max(AC_data.mpl);

ave_frames = cell(num_planes, num_dsets);
std_frames = cell(num_planes, num_dsets);
num_frames = zeros(num_planes, num_dsets);

proc_data = cell(num_dsets,1);
for n_dset = 1:num_dsets
    flist = dir([params.data_dir '\preprocessing\*' AC_data.dset_name{n_dset} '*' AC_data.mouse_tag{n_dset} '*processed_data*']);
    if numel(flist.name)
        proc_data{n_dset} = load([params.data_dir '\preprocessing\' flist.name]);
    end
end

for n_pl = 1:num_planes
    Y = cell(num_dsets,1);
    for n_dset = 1:num_dsets
        flist = dir([params.data_dir '\movies\*' AC_data.dset_name{n_dset} '*' AC_data.mouse_tag{n_dset} '*pl' num2str(n_pl) '*']);
        if numel(flist.name)
            Y_temp = h5read([params.data_dir '\movies\' flist.name], '/mov');
            [d1, d2, ~] = size(Y_temp);
            vid_cuts_trace = logical(proc_data{n_dset}.data.file_cuts_params{n_pl}.vid_cuts_trace);
            Y{n_dset} = zeros(d1, d2, numel(vid_cuts_trace), 'uint16');
            Y{n_dset}(:,:,vid_cuts_trace) = Y_temp;
            
            ave_frames{n_pl, n_dset} = mean(single(Y{n_dset}),3);
            std_frames{n_pl, n_dset} = std(single(Y{n_dset}),0,3);
            num_frames(n_pl, n_dset) = size(Y{n_dset},3);
        end
        
    end
    Y_all = cat(3, Y{:});
    clear Y Y_temp;
    
    
    im1 = cat(2,std_frames{n_pl, :});
    figure;
    imagesc(im1);
    axis equal tight;
    title(sprintf('%s; plane %d; dset 1 - %d', AC_data.mouse_id{n_dset}, n_pl, num_dsets), 'interpreter', 'none');
    
    
    ave_frame_all = mean(single(Y_all),3);
    std_frame_all = std(single(Y_all), 0, 3);
    
    Y_alln = f_s4_remove_mov_bkg(Y_all);
    
    %
    f1 = figure; hold on;
    im1 = imagesc(std_frame_all); axis equal tight;
    im1.Parent.YDir = 'reverse';
    title('click on cell center');
    [n, m] = ginput(1);
    mr = round(m);
    nr = round(n); 
    % idx1[m1, m2, n1, n2]
    idx1 = [(mr-roi_half_size), (mr+roi_half_size); (nr-roi_half_size), (nr+roi_half_size)];
    rectangle('Position', [idx1(2,1), idx1(1,1), abs(diff(idx1(1,:))), abs(diff(idx1(1,:)))], 'EdgeColor', 'r');
    title(sprintf('click to select cell 1'))

    % get mean trace and roi
    roi1 = Y_all(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2),:);
    roi1n = Y_alln(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2),:);
    trace1 = squeeze(mean(mean(roi1,1),2));
    trace1 = trace1 - mean(trace1);
    
    
    [dr1, dr2, Tr] = size(roi1n);
    roi1n2d = reshape(roi1n, dr1*dr2, Tr);

    % get svd version
    %U*S*V'
    [U,S,V] = svd(single(roi1n2d),'econ');
    data_1d = U(:,1)*S(1,1)*V(:,1)';
    spat1 = reshape(U(:,1)*S(1,1)*mean(V(:,1)),dr1, dr2);
    trace1n = mean(U(:,1))*S(1,1)*V(:,1);
    trace1nd = f_smooth_dfdt3(trace1n', 1, 2);
    
    figure;
    subplot(2,2,1);
    imagesc(mean(roi1,3))
    title('mean roi frame')
    subplot(2,2,2);
    imagesc(spat1)
    title('svd extracted roi')
    subplot(2,2,3:4); hold on
    plot(trace1/norm(trace1));
    plot(trace1n/norm(trace1n));
    axis tight;
    title('mean roi trace')
    legend('mean', 'svd extracted')
    sgtitle(sprintf('%s, plane %d', AC_data.mouse_id{1}, n_pl))
    
    
    
    % here need to do stim trig ave
    if sum(strcmpi([AC_data.paradigm], 'behavior'))
        bh_proc = proc_data{strcmpi([AC_data.paradigm], 'behavior')};
        bh_freq_ops = bh_proc.data.stim_params.ops;
        stim_tone_bh = bh_freq_ops.dev_tone_list;
    end
    
    if sum(strcmpi([AC_data.paradigm], 'ammn_stim'))
        ammn_stim_proc = proc_data{strcmpi([AC_data.paradigm], 'ammn_stim')};
        ammn_freq_ops = ammn_stim_proc.data.stim_params.ops;
        stim_run = ammn_freq_ops.stim_trials_volt{1};
        mmn_pat = ammn_freq_ops.MMN_patterns(ammn_freq_ops.paradigm_MMN_pattern(stim_run),:);
        stim_tone_dd = (stim_run - 1)*100 + 70;
        stim_tone_cont = mmn_pat(2);
    end
    
    num_frames_cumsum = cumsum(num_frames(n_pl,:));
    win1 = [5 60];
    
    num_paradigms = numel(AC_data.paradigm);
    for n_pr = 1:num_paradigms
        if sum(strcmpi(AC_data.paradigm{n_pr}, {'behavior', 'ammn'}))
            stim_chan = 'stim type';
        elseif sum(strcmpi(AC_data.paradigm{n_pr}, {'ammn_stim', 'spont_stim'}))
            stim_chan = 'pockel';
        elseif strcmpi(AC_data.paradigm{n_pr}, 'spont')
            stim_chan = '';
        else
            fprintf('paradigm %s not listed\n', AC_data.paradigm{n_pr});
            stim_chan = '';
        end
        
        if numel(stim_chan)
            idx_chan = strcmpi(proc_data{n_pr}.ops.chan_labels, stim_chan);
            stim_times_frame = proc_data{n_pr}.data.stim_times_frame{idx_chan,n_pl};
            
            stim_times_frame(stim_times_frame < win1(1)) = [];
            stim_times_frame((stim_times_frame + win1(2)) > num_frames(n_pr)) = [];
            
            stim_trace2 = trace1n((num_frames_cumsum(n_pr)-num_frames(n_pl,n_pr)+1):num_frames_cumsum(n_pr));
            sorted_trials = squeeze(f_get_stim_trig_resp(stim_trace2', stim_times_frame, win1));
            stim_trace2d = trace1nd((num_frames_cumsum(n_pr)-num_frames(n_pl,n_pr)+1):num_frames_cumsum(n_pr));
            sorted_trialsd = squeeze(f_get_stim_trig_resp(stim_trace2d, stim_times_frame, win1));

            figure; 
            subplot(3,1,1);
            imagesc(sorted_trials'); 
            title('trials all');
            subplot(3,1,2); hold on;
            plot(sorted_trials, 'color', [.6 .6 .6])
            plot(mean(sorted_trials,2), 'Linewidth', 2, 'color', 'k');
            axis tight;
            title('trig ave');
            subplot(3,1,3); hold on;
            plot(sorted_trialsd, 'color', [.6 .6 .6])
            plot(mean(sorted_trialsd,2), 'Linewidth', 2, 'color', 'k');
            axis tight;
            title('trig ave deconvolved');
            sgtitle(stim_chan);
        end
        
        if strcmpi(AC_data.paradigm{n_pr}, 'ammn')
            1
            
        end
        
    end
    
    
    
    idx_stim = strcmpi([AC_data.paradigm], 'ammn_stim');
    if sum(idx_stim)
        
    end
    
    idx_stim = strcmpi([AC_data.paradigm], 'behavior');
    if sum(idx_stim)
        idx_chan = strcmpi(proc_data{idx_stim}.ops.chan_labels, 'stim type');
        stim_times_frame = proc_data{idx_stim}.data.stim_times_frame{idx_chan,n_pl};
        
        stim_trace2 = trace1n((num_frames_cumsum(idx_stim)-num_frames(n_pl,idx_stim)+1):num_frames_cumsum(idx_stim));
        sorted_trials = squeeze(f_get_stim_trig_resp(stim_trace2', stim_times_frame, win1));
        
        figure; 
        subplot(2,1,1);
        imagesc(sorted_trials');
        title('trials all');
        subplot(2,1,2); hold on;
        plot(sorted_trials, 'color', [.6 .6 .6])
        plot(mean(sorted_trials,2), 'Linewidth', 2, 'color', 'k');
        axis tight;
        title('trig ave');
        sgtitle('behavior');
    end
    
    idx_stim = strcmpi([AC_data.paradigm], 'ammn');
    if sum(idx_stim)
        idx_stim = strcmpi(proc_data{idx_stim}.ops.chan_labels, 'stim type');
        stim_times_frame = proc_data{idx_stim}.data.stim_times_frame{idx_chan,n_pl};

    end
    % analyze tuning before + a
    
end

