close all
clear

%%
addpath([pwd '\s1_functions']);
addpath([pwd '\general_functions'])

%%
%ops.file_dir = 'F:\AC_data\caiman_data_dream3\preprocessing';
%ops.file_dir = 'F:\AC_data\caiman_data_echo\preprocessing';

params.data_dir = 'F:\AC_data\extracted_data';

base_onset_win = [0 15];

imprint_stim_key = {'M4463', '11_24_21_pt3', 'spont_stim';...
                    'M4421', '12_24_21b', 170;...
                    'M105', '1_21_22a', 270;...
                    'M4460', '1_2_22b', 270;...
                    'M101', '1_31_22', 170;...
                    'M4463', '11_24_21_pt4', 170;...
                    'M4465', 'a_2_22a', 270;...
                    'M124', '4_11_22', 5;...
                    'M125', '4_26_22', 170;...
                    'M108', '2_4_22a', 3;...
                    'M142', '6_11_22b', 270;...
                    'M166', '6_20_22', 170;
                    'M166', '6_20_22_pt2', 170};
                

%%
flist = dir([params.data_dir '\*_manual_roi.mat']);

%% load all data
num_fl = numel(flist);
fl_AC_data = cell(num_fl, 1);
fl_proc_data = cell(num_fl, 1);
fl_cell_traces = cell(num_fl, 1);
fl_cell_type_idx = cell(num_fl, 1);

for n_fl = 1:num_fl
    load_data = load([params.data_dir '\' flist(n_fl).name]);
    fl_AC_data{n_fl} = load_data.data.AC_data;
    fl_proc_data{n_fl} = load_data.data.proc_data;
    fl_cell_traces{n_fl} = load_data.data.cell_ca;
    fl_cell_type_idx{n_fl} = load_data.data.cell_run_id;
end

fl_cell_traces_prs = cell(num_fl, 1);
fl_cell_traces_prsd = cell(num_fl, 1);
fl_cell_type_idx_cat = cell(num_fl, 1);
for n_fl = 1:num_fl
    fl_cell_type_idx_cat{n_fl} = cat(1,fl_cell_type_idx{n_fl}{:});
    
    temp_trace = cat(1,fl_cell_traces{n_fl}{:});
    temp_trace_d = f_smooth_dfdt3(temp_trace,1, 2, 1);
    
    num_dsets = numel(fl_AC_data{n_fl}.paradigm);
    fl_cell_traces_prs{n_fl} = cell(num_dsets,1);
    fl_cell_traces_prsd{n_fl} = cell(num_dsets,1);
    
    cur_start = 1;
    for n_dset = 1:num_dsets
        num_vol = fl_proc_data{n_fl}{n_dset}.data.frame_data.num_volumes_linear;
        cur_end = cur_start +  num_vol - 1;
        fl_cell_traces_prsd{n_fl}{n_dset} = temp_trace_d(:,cur_start:cur_end);
        fl_cell_traces_prs{n_fl}{n_dset} = temp_trace(:,cur_start:cur_end);
        cur_start = cur_end + 1;
    end
    
end
fl_cell_type_stim = fl_cell_type_idx_cat;

%% get stim resp cells

bh_stim = false(num_fl,1);
ammn_stim = false(num_fl,1);
spont_stim = false(num_fl,1);
for n_fl = 1:num_fl
    AC_data = fl_AC_data{n_fl};
    if sum(strcmpi(AC_data.paradigm, 'ammn_stim'))
        ammn_stim(n_fl) = 1;
    end
    if sum(strcmpi(AC_data.paradigm, 'spont_stim'))
        spont_stim(n_fl) = 1;
    end
    if sum(strcmpi(AC_data.mouse_tag, '11_24_21_pt3'))
        spont_stim(n_fl) = 1;
    end
    if sum(strcmpi(AC_data.paradigm, 'behavior'))
        bh_stim(n_fl) = 1;
    end
    
end
ammn_stim_nobh = ammn_stim;
ammn_stim_nobh(bh_stim) = 0;
%% look throguh all targeted cells and throw away that did not activate

z_thresh = 4;
peak_bin_size = 3;
num_samp = 1000;
peak_prcntle = normcdf(z_thresh)*100;


trials_analysis = [1:10, 170, 270];
trials_sample = [1:10, 170, 270];

num_tt = numel(trials_analysis);

win2 = [5 12];
num_t = sum(win2);

fl_resp_cells = cell(num_fl,1);
fl_resp_thresh = cell(num_fl,1);
fl_resp_cells_z = cell(num_fl,1);
for n_fl = 1:num_fl
    AC_data = fl_AC_data{n_fl};
    
    num_parad = numel(fl_cell_traces_prsd{n_fl});
    
    fl_resp_cells{n_fl} = cell(num_parad,1);
    fl_resp_thresh{n_fl} = cell(num_parad,1);
    fl_resp_cells_z{n_fl} = cell(num_parad,1);
    for n_pr = 1:num_parad
        if sum(strcmpi(AC_data.paradigm{n_pr}, {'ammn', 'ammn_stim'}))
            % compute tuning
            firing_rate = fl_cell_traces_prsd{n_fl}{n_pr};
            proc_data2 = fl_proc_data{n_fl}{n_pr};
            
            num_cells = size(firing_rate,1);
            
            stim_times = proc_data2.data.stim_times_frame{1,1};
            trial_types = proc_data2.data.trial_types;
            
            
            trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, win2);
            
            num_trial_per_stim = min(sum(trial_types == trials_sample,1));
            
            
            stat_idx = logical(sum(trial_types == trials_sample,2));
            pop_stim_times = stim_times(stat_idx);
            trial_data_sort_stat = trial_data_sort(:,:,stat_idx);
            trial_data_stat_mean = mean(trial_data_sort_stat,3);
            trial_data_stat_sem = std(trial_data_sort_stat,[],3)/sqrt(num_trial_per_stim-1);

            num_trials = numel(pop_stim_times);
            
            peak_vals = nan(num_cells, num_tt);
            peak_locs = nan(num_cells, num_tt);
            for n_tt = 1:num_tt
                trial_data_sort2 = trial_data_sort(:,:, trial_types==trials_sample(n_tt));
                if ~isempty(trial_data_sort2)
                    [peak_vals(:,n_tt), peak_locs(:,n_tt)] = f_get_trial_peak(mean(trial_data_sort2,3), peak_bin_size);
                end
            end
            
            if 1
                samp_peak_vals = zeros(num_cells, num_samp);
                samp_peak_locs = zeros(num_cells, num_samp);
                for n_cell = 1:num_cells
                    samp_idx = randsample(num_trials, num_trial_per_stim*num_samp, 1);
                    samp_trial_data_sort = mean(reshape(trial_data_sort(n_cell, :, samp_idx), [num_t, num_samp, num_trial_per_stim]),3)';
                    [samp_peak_vals(n_cell,:), samp_peak_locs(n_cell,:)] = f_get_trial_peak(samp_trial_data_sort, peak_bin_size);
                end
                
                resp_thresh = repmat(prctile(samp_peak_vals', peak_prcntle)', [1 num_t]);
                resp_cells = peak_vals>resp_thresh(:,1);
                resp_cells_z = peak_vals./resp_thresh(:,1);
                
            else
                resp_thresh = zeros(num_cells, num_t);
                resp_cells = zeros(num_cells, num_tt);
                for n_cell = 1:num_cells
                    trial_data1 = squeeze(trial_data_sort_stat(n_cell,:,:));
                    resp_thresh(n_cell,:) = mean(trial_data1,2) + z_thresh*std(trial_data1,[],2)/sqrt(num_trial_per_stim-1);
                    for n_tt = 1:num_tt
                        resp_cells(n_cell, n_tt) = peak_vals(n_cell, n_tt) > resp_thresh(n_cell,peak_locs(n_cell,n_tt));
                    end
                end
            end
            fl_resp_cells{n_fl}{n_pr} = resp_cells;
            fl_resp_thresh{n_fl}{n_pr} = resp_thresh;
            fl_resp_cells_z{n_fl}{n_pr} = resp_cells_z;
        end
    end
end

%% find list of tuned cells in ammn
win1 = [10 10];
for n_fl = 1:num_fl
    AC_data = fl_AC_data{n_fl};
    
    stim_idx = find(strcmpi([AC_data.paradigm], {'ammn_stim'}),1);
    stim_chan = 'pockel';
    fl_cell_type_stim{n_fl} = false(numel(fl_cell_type_idx_cat{n_fl}),1);
    if sum(stim_idx)
    
        proc_data2 = fl_proc_data{n_fl}{stim_idx};
        
        idx_chan = strcmpi(proc_data2.ops.chan_labels, stim_chan);
        stim_times_frame = proc_data2.data.stim_times_frame{idx_chan,1};
        
        targ_cells= find(fl_cell_type_idx_cat{n_fl}==1);
        
        trace1n = fl_cell_traces_prsd{n_fl}{stim_idx}(fl_cell_type_idx_cat{n_fl}==1,:);
        trace1nd = fl_cell_traces_prsd{n_fl}{stim_idx}(fl_cell_type_idx_cat{n_fl}==1,:);
        
        [num_cells, num_T] = size(trace1n);
        
        
        stim_times_frame(stim_times_frame < win1(1)) = [];
        stim_times_frame((stim_times_frame + win1(2)) > num_T) = [];

        
        %stim_trace2 = trace1n((num_frames_cumsum(n_pr)-num_vol_all(n_pr)+1):num_frames_cumsum(n_pr));
        sorted_trials = f_get_stim_trig_resp(trace1n, stim_times_frame, win1);
        %tim_trace2d = trace1nd((num_frames_cumsum(n_pr)-num_vol_all(n_pr)+1):num_frames_cumsum(n_pr));
        sorted_trialsd = f_get_stim_trig_resp(trace1nd, stim_times_frame, win1);
        
        tr_mean_sorted_trialsd = mean(sorted_trialsd,3);
        
        base1 = mean(tr_mean_sorted_trialsd(:,1:win1(1)),2);
        resp1 = mean(tr_mean_sorted_trialsd(:,(win1(1)+1):end),2);
        stim_tuning = (resp1./base1);
        
        
        fl_cell_type_stim{n_fl}(targ_cells(stim_tuning > 2)) = 1;
        
        if 0
            for n_cell = 1:num_cells
                if stim_tuning(n_cell) > 1.3
                    cell_sorted_trials = squeeze(sorted_trials(n_cell,:,:));
                    cell_sorted_trialsd = squeeze(sorted_trialsd(n_cell,:,:));
                    figure; 
                    subplot(3,1,1);
                    imagesc(cell_sorted_trials'); 
                    title('trials all');
                    subplot(3,1,2); hold on;
                    plot(cell_sorted_trials, 'color', [.6 .6 .6])
                    plot(mean(cell_sorted_trials,2), 'Linewidth', 2, 'color', 'k');
                    axis tight;
                    title('trig ave');
                    subplot(3,1,3); hold on;
                    plot(cell_sorted_trialsd, 'color', [.6 .6 .6])
                    plot(mean(cell_sorted_trialsd,2), 'Linewidth', 2, 'color', 'k');
                    axis tight;
                    title('trig ave deconvolved');
                    sgtitle(sprintf('%s; cell %d, %s; opto tuning %.2f',stim_chan, n_cell, AC_data.paradigm{stat_idx}, stim_tuning(n_cell)));
                end
            end
        end
    end
end



%% first ammn no bh

fl_cell_traces_prsd2 = fl_cell_traces_prsd(ammn_stim_nobh);
fl_proc_data2 = fl_proc_data(ammn_stim_nobh);
fl_AC_data2 = fl_AC_data(ammn_stim_nobh);
fl_cell_type_stim2 = fl_cell_type_stim(ammn_stim_nobh);
fl_cell_type_idx_cat2 = fl_cell_type_idx_cat(ammn_stim_nobh);
fl_resp_cells2 = fl_resp_cells(ammn_stim_nobh);

num_fl = sum(ammn_stim_nobh);

params.base_onset_win = [0 15];

params.corr_method = 'cosine';
params.plot_sing_tr = 0;
params.plot_means = 0;
params.plot_corr = 1;
params.plot_rast = 0;
params.plot_sig_corr = 1;
params.norm_to_all_tr = 1;

params.set_mean_max = 0.7;
params.set_corr_ylim = [0 .5];
params.set_sig_corr_ylim = [0 .5];


corr_sig_ammn_nobh = cell(num_fl,1);
corr_noise_ammn_nobh = cell(num_fl,1);
corr_sig_ammn_nobh_samedd = cell(num_fl,1);
corr_noise_ammn_nobh_samedd = cell(num_fl,1);
corr_sig_ammn_nobh_otherdd = cell(num_fl,1);
corr_noise_ammn_nobh_otherdd = cell(num_fl,1);
rest_corr_all = cell(num_fl, 1);
rest_corr_to_other = cell(num_fl, 1);
for n_fl = 1:num_fl
    if sum(fl_cell_type_stim2{n_fl})
        
        ammn_idx = find(strcmpi([fl_AC_data2{n_fl}.paradigm], {'ammn'}),1);
        
        stat_idx = strcmpi([fl_AC_data2{n_fl}.paradigm], {'ammn_stim'});
        if sum(stat_idx) > 1
            idx2 = find(stat_idx,1);
            fl_cell_traces_prsd2{n_fl}(idx2) = [];
            fl_proc_data2{n_fl}(idx2) = [];
            fl_AC_data2{n_fl}(idx2,:) = [];
        end
        
        disp(fl_AC_data2{n_fl}.mouse_id{1})
        params.paradigm_to_check = {'ammn', 'ammn_stim'};
        params.cell_type_idx = fl_cell_type_stim2{n_fl};
        
        idx1 = strcmpi(fl_AC_data2{n_fl}.mouse_tag{1}, imprint_stim_key(:,2));
        params.trials_to_use_for_norm = [1:10 170 270];
         %                     170, 270,...
        %                      4, 7,... %,...
        %                      201:206};
                             %1, 2, 3, 5, 6, 8, 9,...
                             %1 2 3 4 5 6 8 9 10,...
                             %201:206};%,...
                             %101:106};
                             
        params.corr_trials_check = {imprint_stim_key{idx1,3}};
        
        params_samedd = params;
        params_samedd.cell_type_idx = fl_resp_cells2{n_fl}{ammn_idx}(:,trials_analysis==params_samedd.corr_trials_check{1});
        
        params_otherdd = params;
        params_otherdd.corr_trials_check = {440 - imprint_stim_key{idx1,3}};
        params_otherdd.cell_type_idx = fl_resp_cells2{n_fl}{ammn_idx}(:,trials_analysis==params_otherdd.corr_trials_check{1});
        
        data_out = f_s5_sig_noise_corr_within(fl_cell_traces_prsd2{n_fl}, fl_proc_data2{n_fl}, fl_AC_data2{n_fl}, params);
        corr_sig_ammn_nobh{n_fl} = data_out.sig_corr.pairs_all;
        corr_noise_ammn_nobh{n_fl} = data_out.noise_corr.pairs_all;
        sig_corr_labels = data_out.sig_corr.labels_ammn;
        noise_corr_labels = data_out.noise_corr.labels_ammn;
        
        if sum(params_samedd.cell_type_idx)>1
            data_out = f_s5_sig_noise_corr_within(fl_cell_traces_prsd2{n_fl}, fl_proc_data2{n_fl}, fl_AC_data2{n_fl}, params_samedd);
            corr_sig_ammn_nobh_samedd{n_fl} = data_out.sig_corr.pairs_all;
            corr_noise_ammn_nobh_samedd{n_fl} = data_out.noise_corr.pairs_all;
        end
        if sum(params_otherdd.cell_type_idx)>1
            data_out = f_s5_sig_noise_corr_within(fl_cell_traces_prsd2{n_fl}, fl_proc_data2{n_fl}, fl_AC_data2{n_fl}, params_otherdd);
            corr_sig_ammn_nobh_otherdd{n_fl} = data_out.sig_corr.pairs_all;
            corr_noise_ammn_nobh_otherdd{n_fl} = data_out.noise_corr.pairs_all;
        end
        num_cells = sum(fl_cell_type_stim2{n_fl});
        rest_files = find(strcmpi([fl_AC_data2{n_fl}.paradigm], 'spont'));

        pairs = nchoosek(1:num_cells,2);

        
        SI_all = cell(numel(rest_files), 1);

        for n_file = 1:numel(rest_files)
            SI_all{n_file} = 1-f_pdist_YS(fl_cell_traces_prsd2{n_fl}{rest_files(n_file)}(fl_cell_type_stim2{n_fl},:), 'cosine');
        end

        SI_all_cat = cat(3,SI_all{:});
        paris_all = zeros(size(pairs,1), numel(rest_files));
        for n_pair = 1:size(pairs,1)
            paris_all(n_pair, :) =  squeeze(SI_all_cat(pairs(n_pair,1),pairs(n_pair,2),:));
        end
        rest_corr_all{n_fl} = paris_all;
        
        % noe mabe imprinted to all other 
        
        SI_all = cell(numel(rest_files), 1);

        for n_file = 1:numel(rest_files)
            vec1 = fl_cell_traces_prsd2{n_fl}{rest_files(n_file)}(fl_cell_type_stim2{n_fl},:);
            vec2 = fl_cell_traces_prsd2{n_fl}{rest_files(n_file)}(fl_cell_type_idx_cat2{n_fl}==3,:);
            tempsi = 1-pdist2(vec1,vec2, 'cosine');
            SI_all{n_file} = tempsi(:);
        end
        rest_corr_to_other{n_fl} = cat(2,SI_all{:});

        
%         corr_all_cat = cat(3,corr_all{:});
%         paris_all = zeros(size(pairs,1), numel(rest_files));
%         for n_pair = 1:size(pairs,1)
%             paris_all(n_pair, :) =  squeeze(corr_all_cat(pairs(n_pair,1),pairs(n_pair,2),rest_files));
%         end
%         [~,p] = ttest(paris_all(:,1),paris_all(:,end));
%         figure; hold on;
%         plot(labels_rest, paris_all, 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
%         errorbar(labels_rest, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
%         ylabel('pairwise correlations')
%         title(sprintf('rest pre-post pval=%.2f', p));
% 

        
        
    end
end

%% plot nobh ammn corr
corr_sig_ammn_nobh2 = cat(1,corr_sig_ammn_nobh{:});
corr_noise_ammn_nobh2 = cat(1,corr_noise_ammn_nobh{:});
rest_corr_all2  = cat(1,rest_corr_all{:});
rest_corr_to_other2 = cat(1,rest_corr_to_other{:});

[~,p] = ttest(corr_sig_ammn_nobh2(:,1),corr_sig_ammn_nobh2(:,end));
figure; hold on;
plot(sig_corr_labels, corr_sig_ammn_nobh2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh2,1),std(corr_sig_ammn_nobh2,[],1)./sqrt(size(corr_sig_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('signal correlations ammn pre-post; pval=%.2f', p));


[~,p] = ttest(corr_noise_ammn_nobh2(:,1),corr_noise_ammn_nobh2(:,end));
figure; hold on;
plot(noise_corr_labels, corr_noise_ammn_nobh2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(noise_corr_labels, mean(corr_noise_ammn_nobh2,1),std(corr_noise_ammn_nobh2,[],1)./sqrt(size(corr_noise_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('noise correlations ammn pre-post; pval=%.2f', p));

labels_cell_rest = {'pre', 'post'};
labels_rest = reordercats(categorical(labels_cell_rest),labels_cell_rest);


[~,p] = ttest(rest_corr_all2(:,1),rest_corr_all2(:,end));
figure; hold on;
plot(labels_rest, rest_corr_all2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(rest_corr_all2,1),std(rest_corr_all2,[],1)./sqrt(size(rest_corr_all2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('rest pre-post pval=%.2f', p));

rest_corr_to_other3 = rest_corr_to_other2(~isnan(rest_corr_to_other2(:,2)),:);
[~,p] = ttest(rest_corr_to_other3(:,1),rest_corr_to_other3(:,end));
figure; hold on;
plot(labels_rest, rest_corr_to_other3', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(rest_corr_to_other3,1),std(rest_corr_to_other3,[],1)./sqrt(size(rest_corr_to_other3,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('rest ot other pre-post pval=%.2f', p));
ylim([.2 .35])


rest_corr_to_other{n_fl} = SI_all;


%%

corr_sig_ammn_nobh_samedd2 = cat(1,corr_sig_ammn_nobh_samedd{:});
corr_noise_ammn_nobh_samedd2 = cat(1,corr_noise_ammn_nobh_samedd{:});

corr_sig_ammn_nobh_otherdd2 = cat(1,corr_sig_ammn_nobh_otherdd{:});
corr_noise_ammn_nobh_otherdd2 = cat(1,corr_noise_ammn_nobh_otherdd{:});



[~,p] = ttest(corr_sig_ammn_nobh_samedd2(:,1),corr_sig_ammn_nobh_samedd2(:,end));
figure; hold on;
plot(sig_corr_labels, corr_sig_ammn_nobh_samedd2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh_samedd2,1),std(corr_sig_ammn_nobh_samedd2,[],1)./sqrt(size(corr_sig_ammn_nobh_samedd2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('same dd signal correlations ammn pre-post; pval=%.2f', p));

[~,p] = ttest(corr_sig_ammn_nobh_otherdd2(:,1),corr_sig_ammn_nobh_otherdd2(:,end));
figure; hold on;
plot(sig_corr_labels, corr_sig_ammn_nobh_otherdd2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh_otherdd2,1),std(corr_sig_ammn_nobh_otherdd2,[],1)./sqrt(size(corr_sig_ammn_nobh_otherdd2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('other dd signal correlations ammn pre-post; pval=%.2f', p));



figure; hold on;
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh2,1),std(corr_sig_ammn_nobh2,[],1)./sqrt(size(corr_sig_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2)
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh_samedd2,1),std(corr_sig_ammn_nobh_samedd2,[],1)./sqrt(size(corr_sig_ammn_nobh_samedd2,1)-1), 'o-', 'Linewidth', 2)
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh_otherdd2,1),std(corr_sig_ammn_nobh_otherdd2,[],1)./sqrt(size(corr_sig_ammn_nobh_otherdd2,1)-1), 'o-', 'Linewidth', 2)
legend('stim cells', 'same dd', 'other dd')
title('signal corr')

figure; hold on;
errorbar(sig_corr_labels, mean(corr_noise_ammn_nobh2,1),std(corr_noise_ammn_nobh2,[],1)./sqrt(size(corr_noise_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2)
errorbar(sig_corr_labels, mean(corr_noise_ammn_nobh_samedd2,1),std(corr_noise_ammn_nobh_samedd2,[],1)./sqrt(size(corr_noise_ammn_nobh_samedd2,1)-1), 'o-', 'Linewidth', 2)
errorbar(sig_corr_labels, mean(corr_noise_ammn_nobh_otherdd2,1),std(corr_noise_ammn_nobh_otherdd2,[],1)./sqrt(size(corr_noise_ammn_nobh_otherdd2,1)-1), 'o-', 'Linewidth', 2)
legend('stim cells', 'same dd', 'other dd')
title('noise corr')




%% now with bh

fl_cell_traces_prsd2 = fl_cell_traces_prsd(bh_stim);
fl_proc_data2 = fl_proc_data(bh_stim);
fl_AC_data2 = fl_AC_data(bh_stim);
fl_cell_type_stim2 = fl_cell_type_stim(bh_stim);

num_fl = sum(bh_stim);

params.base_onset_win = [0 15];

params.corr_method = 'cosine';
params.plot_sing_tr = 0;
params.plot_means = 0;
params.plot_corr = 1;
params.plot_rast = 0;
params.plot_sig_corr = 1;
params.norm_to_all_tr = 1;

params.set_mean_max = 0.7;
params.set_corr_ylim = [0 .5];
params.set_sig_corr_ylim = [0 .5];


corr_sig_ammn_nobh = cell(num_fl,1);
corr_noise_ammn_nobh = cell(num_fl,1);
rest_corr_all = cell(num_fl, 1);

for n_fl = 1:num_fl
    if sum(fl_cell_type_stim2{n_fl})
        
        stat_idx = strcmpi([fl_AC_data2{n_fl}.paradigm], {'ammn_stim'});
        if sum(stat_idx) > 1
            idx2 = find(stat_idx,1);
            fl_cell_traces_prsd2{n_fl}(idx2) = [];
            fl_proc_data2{n_fl}(idx2) = [];
            fl_AC_data2{n_fl}(idx2,:) = [];
        end
        
        disp(fl_AC_data2{n_fl}.mouse_id{1})
        params.paradigm_to_check = {'ammn', 'ammn_stim'};
        params.cell_type_idx = fl_cell_type_stim2{n_fl};
        
        idx1 = strcmpi(fl_AC_data2{n_fl}.mouse_tag{1}, imprint_stim_key(:,2));
        
        params.corr_trials_check = {imprint_stim_key{idx1,3}};
        params.trials_to_use_for_norm = [1:10 170 270];
        %                     170, 270,...
        %                      4, 7,... %,...
        %                      201:206};
                             %1, 2, 3, 5, 6, 8, 9,...
                             %1 2 3 4 5 6 8 9 10,...
                             %201:206};%,...
                             %101:106};


        data_out = f_s5_sig_noise_corr_within(fl_cell_traces_prsd2{n_fl}, fl_proc_data2{n_fl}, fl_AC_data2{n_fl}, params);
        corr_sig_ammn_nobh{n_fl} = data_out.sig_corr.pairs_all;
        corr_noise_ammn_nobh{n_fl} = data_out.noise_corr.pairs_all;
        sig_corr_labels = data_out.sig_corr.labels_ammn;
        noise_corr_labels = data_out.noise_corr.labels_ammn;
        
        
        num_cells = sum(fl_cell_type_stim2{n_fl});
        rest_files = find(strcmpi([fl_AC_data2{n_fl}.paradigm], 'spont'));

        pairs = nchoosek(1:num_cells,2);

        
        SI_all = cell(numel(rest_files), 1);

        for n_file = 1:numel(rest_files)
            SI_all{n_file} = 1-f_pdist_YS(fl_cell_traces_prsd2{n_fl}{rest_files(n_file)}(fl_cell_type_stim2{n_fl},:), 'cosine');
        end

        SI_all_cat = cat(3,SI_all{:});
        paris_all = zeros(size(pairs,1), numel(rest_files));
        for n_pair = 1:size(pairs,1)
            paris_all(n_pair, :) =  squeeze(SI_all_cat(pairs(n_pair,1),pairs(n_pair,2),:));
        end
        rest_corr_all{n_fl} = paris_all;
        
        
%         corr_all_cat = cat(3,corr_all{:});
%         paris_all = zeros(size(pairs,1), numel(rest_files));
%         for n_pair = 1:size(pairs,1)
%             paris_all(n_pair, :) =  squeeze(corr_all_cat(pairs(n_pair,1),pairs(n_pair,2),rest_files));
%         end
%         [~,p] = ttest(paris_all(:,1),paris_all(:,end));
%         figure; hold on;
%         plot(labels_rest, paris_all, 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
%         errorbar(labels_rest, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
%         ylabel('pairwise correlations')
%         title(sprintf('rest pre-post pval=%.2f', p));
% 

        
        
    end
end

%% plot wit hbh
corr_sig_ammn_nobh2 = cat(1,corr_sig_ammn_nobh{:});
corr_noise_ammn_nobh2 = cat(1,corr_noise_ammn_nobh{:});
rest_corr_all2  = cat(1,rest_corr_all{:});

[~,p] = ttest(corr_sig_ammn_nobh2(:,1),corr_sig_ammn_nobh2(:,end));
figure; hold on;
plot(sig_corr_labels, corr_sig_ammn_nobh2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh2,1),std(corr_sig_ammn_nobh2,[],1)./sqrt(size(corr_sig_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('signal correlations ammn pre-post; pval=%.2f', p));


[~,p] = ttest(corr_noise_ammn_nobh2(:,1),corr_noise_ammn_nobh2(:,end));
figure; hold on;
plot(noise_corr_labels, corr_noise_ammn_nobh2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(noise_corr_labels, mean(corr_noise_ammn_nobh2,1),std(corr_noise_ammn_nobh2,[],1)./sqrt(size(corr_noise_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('noise correlations ammn pre-post; pval=%.2f', p));

labels_cell_rest = {'pre', 'post'};
labels_rest = reordercats(categorical(labels_cell_rest),labels_cell_rest);


[~,p] = ttest(rest_corr_all2(:,1),rest_corr_all2(:,end));
figure; hold on;
plot(labels_rest, rest_corr_all2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(rest_corr_all2,1),std(rest_corr_all2,[],1)./sqrt(size(rest_corr_all2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('rest pre-post pval=%.2f', p));

%%


fl_cell_traces_prsd2 = fl_cell_traces_prsd(bh_stim);
fl_proc_data2 = fl_proc_data(bh_stim);
fl_AC_data2 = fl_AC_data(bh_stim);
fl_cell_type_stim2 = fl_cell_type_stim(bh_stim);

num_fl = sum(bh_stim);

params.base_onset_win = [0 15];

params.corr_method = 'cosine';
params.plot_sing_tr = 0;
params.plot_means = 0;
params.plot_corr = 1;
params.plot_rast = 0;
params.plot_sig_corr = 1;
params.norm_to_all_tr = 1;

params.set_mean_max = 0.7;
params.set_corr_ylim = [0 .5];
params.set_sig_corr_ylim = [0 .5];


corr_sig_ammn_nobh = cell(num_fl,1);
corr_noise_ammn_nobh = cell(num_fl,1);
rest_corr_all = cell(num_fl, 1);

for n_fl = 1:num_fl
    if sum(fl_cell_type_stim2{n_fl})
        
    end
end

%% now rest stim
fl_cell_traces_prsd2 = fl_cell_traces_prsd(spont_stim);
fl_proc_data2 = fl_proc_data(spont_stim);
fl_AC_data2 = fl_AC_data(spont_stim);
fl_cell_type_stim2 = fl_cell_type_stim(spont_stim);

num_fl = sum(spont_stim);

params.base_onset_win = [0 15];

params.corr_method = 'cosine';
params.plot_sing_tr = 0;
params.plot_means = 0;
params.plot_corr = 1;
params.plot_rast = 0;
params.plot_sig_corr = 1;
params.norm_to_all_tr = 1;

params.set_mean_max = 0.7;
params.set_corr_ylim = [0 .5];
params.set_sig_corr_ylim = [0 .5];


corr_sig_ammn_nobh = cell(num_fl,1);
corr_noise_ammn_nobh = cell(num_fl,1);
rest_corr_all = cell(num_fl, 1);

for n_fl = 1:num_fl
    if sum(fl_cell_type_stim2{n_fl})
        
        stat_idx = strcmpi([fl_AC_data2{n_fl}.paradigm], {'ammn_stim'});
        if sum(stat_idx) > 1
            idx2 = find(stat_idx,1);
            fl_cell_traces_prsd2{n_fl}(idx2) = [];
            fl_proc_data2{n_fl}(idx2) = [];
            fl_AC_data2{n_fl}(idx2,:) = [];
        end
        
        disp(fl_AC_data2{n_fl}.mouse_id{1})
        params.paradigm_to_check = {'ammn', 'ammn_stim'};
        params.cell_type_idx = fl_cell_type_stim2{n_fl};
        
        idx1 = strcmpi(fl_AC_data2{n_fl}.mouse_tag{1}, imprint_stim_key(:,2));
        
        params.corr_trials_check = {imprint_stim_key{idx1,3}};
        params.trials_to_use_for_norm = [1:10 170 270];
        %                     170, 270,...
        %                      4, 7,... %,...
        %                      201:206};
                             %1, 2, 3, 5, 6, 8, 9,...
                             %1 2 3 4 5 6 8 9 10,...
                             %201:206};%,...
                             %101:106};


        data_out = f_s5_sig_noise_corr_within(fl_cell_traces_prsd2{n_fl}, fl_proc_data2{n_fl}, fl_AC_data2{n_fl}, params);
        corr_sig_ammn_nobh{n_fl} = data_out.sig_corr.pairs_all;
        corr_noise_ammn_nobh{n_fl} = data_out.noise_corr.pairs_all;
        sig_corr_labels = data_out.sig_corr.labels_ammn;
        noise_corr_labels = data_out.noise_corr.labels_ammn;
        
        
        num_cells = sum(fl_cell_type_stim2{n_fl});
        rest_files = find(strcmpi([fl_AC_data2{n_fl}.paradigm], 'spont'));

        pairs = nchoosek(1:num_cells,2);

        
        SI_all = cell(numel(rest_files), 1);

        for n_file = 1:numel(rest_files)
            SI_all{n_file} = 1-f_pdist_YS(fl_cell_traces_prsd2{n_fl}{rest_files(n_file)}(fl_cell_type_stim2{n_fl},:), 'cosine');
        end

        SI_all_cat = cat(3,SI_all{:});
        paris_all = zeros(size(pairs,1), numel(rest_files));
        for n_pair = 1:size(pairs,1)
            paris_all(n_pair, :) =  squeeze(SI_all_cat(pairs(n_pair,1),pairs(n_pair,2),:));
        end
        rest_corr_all{n_fl} = paris_all;
        
        
%         corr_all_cat = cat(3,corr_all{:});
%         paris_all = zeros(size(pairs,1), numel(rest_files));
%         for n_pair = 1:size(pairs,1)
%             paris_all(n_pair, :) =  squeeze(corr_all_cat(pairs(n_pair,1),pairs(n_pair,2),rest_files));
%         end
%         [~,p] = ttest(paris_all(:,1),paris_all(:,end));
%         figure; hold on;
%         plot(labels_rest, paris_all, 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
%         errorbar(labels_rest, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
%         ylabel('pairwise correlations')
%         title(sprintf('rest pre-post pval=%.2f', p));
% 

        
        
    end
end

%% plot now rest stim
corr_sig_ammn_nobh2 = cat(1,corr_sig_ammn_nobh{:});
corr_noise_ammn_nobh2 = cat(1,corr_noise_ammn_nobh{:});
rest_corr_all2  = cat(1,rest_corr_all{:});

[~,p] = ttest(corr_sig_ammn_nobh2(:,1),corr_sig_ammn_nobh2(:,end));
figure; hold on;
plot(sig_corr_labels, corr_sig_ammn_nobh2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh2,1),std(corr_sig_ammn_nobh2,[],1)./sqrt(size(corr_sig_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('signal correlations ammn pre-post; pval=%.2f', p));


[~,p] = ttest(corr_noise_ammn_nobh2(:,1),corr_noise_ammn_nobh2(:,end));
figure; hold on;
plot(noise_corr_labels, corr_noise_ammn_nobh2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(noise_corr_labels, mean(corr_noise_ammn_nobh2,1),std(corr_noise_ammn_nobh2,[],1)./sqrt(size(corr_noise_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('noise correlations ammn pre-post; pval=%.2f', p));

labels_cell_rest = {'pre', 'post'};
labels_rest = reordercats(categorical(labels_cell_rest),labels_cell_rest);


[~,p] = ttest(rest_corr_all2(:,1),rest_corr_all2(:,end));
figure; hold on;
plot(labels_rest, rest_corr_all2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(rest_corr_all2,1),std(rest_corr_all2,[],1)./sqrt(size(rest_corr_all2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('rest pre-post pval=%.2f', p));


%%
fl_cell_traces_prsd2 = fl_cell_traces_prsd(ammn_stim_nobh);
fl_proc_data2 = fl_proc_data(ammn_stim_nobh);
fl_AC_data2 = fl_AC_data(ammn_stim_nobh);
fl_cell_type_stim2 = fl_cell_type_stim(ammn_stim_nobh);

num_fl = sum(ammn_stim_nobh);

params.base_onset_win = [0 15];

params.corr_method = 'cosine';
params.plot_sing_tr = 0;
params.plot_means = 0;
params.plot_corr = 1;
params.plot_rast = 0;
params.plot_sig_corr = 1;
params.norm_to_all_tr = 1;

params.set_mean_max = 0.7;
params.set_corr_ylim = [0 .5];
params.set_sig_corr_ylim = [0 .5];


corr_sig_ammn_nobh = cell(num_fl,1);
corr_noise_ammn_nobh = cell(num_fl,1);
rest_corr_all = cell(num_fl, 1);

for n_fl = 1:num_fl
    if sum(fl_cell_type_stim2{n_fl})
        
        stat_idx = strcmpi([fl_AC_data2{n_fl}.paradigm], {'ammn_stim'});
        if sum(stat_idx) > 1
            idx2 = find(stat_idx,1);
            fl_cell_traces_prsd2{n_fl}(idx2) = [];
            fl_proc_data2{n_fl}(idx2) = [];
            fl_AC_data2{n_fl}(idx2,:) = [];
        end
        
        disp(fl_AC_data2{n_fl}.mouse_id{1})
        params.paradigm_to_check = {'ammn', 'ammn_stim'};
        params.cell_type_idx = fl_cell_type_stim2{n_fl};
        
        idx1 = strcmpi(fl_AC_data2{n_fl}.mouse_tag{1}, imprint_stim_key(:,2));
        
        params.corr_trials_check = {imprint_stim_key{idx1,3}};
        params.trials_to_use_for_norm = [1:10 170 270];
        %                     170, 270,...
        %                      4, 7,... %,...
        %                      201:206};
                             %1, 2, 3, 5, 6, 8, 9,...
                             %1 2 3 4 5 6 8 9 10,...
                             %201:206};%,...
                             %101:106};


        data_out = f_s5_sig_noise_corr_within(fl_cell_traces_prsd2{n_fl}, fl_proc_data2{n_fl}, fl_AC_data2{n_fl}, params);
        corr_sig_ammn_nobh{n_fl} = data_out.sig_corr.pairs_all;
        corr_noise_ammn_nobh{n_fl} = data_out.noise_corr.pairs_all;
        sig_corr_labels = data_out.sig_corr.labels_ammn;
        noise_corr_labels = data_out.noise_corr.labels_ammn;
        
        
        num_cells = sum(fl_cell_type_stim2{n_fl});
        rest_files = find(strcmpi([fl_AC_data2{n_fl}.paradigm], 'spont'));

        pairs = nchoosek(1:num_cells,2);

        
        SI_all = cell(numel(rest_files), 1);

        for n_file = 1:numel(rest_files)
            SI_all{n_file} = 1-f_pdist_YS(fl_cell_traces_prsd2{n_fl}{rest_files(n_file)}(fl_cell_type_stim2{n_fl},:), 'cosine');
        end

        SI_all_cat = cat(3,SI_all{:});
        paris_all = zeros(size(pairs,1), numel(rest_files));
        for n_pair = 1:size(pairs,1)
            paris_all(n_pair, :) =  squeeze(SI_all_cat(pairs(n_pair,1),pairs(n_pair,2),:));
        end
        rest_corr_all{n_fl} = paris_all;
        
        
%         corr_all_cat = cat(3,corr_all{:});
%         paris_all = zeros(size(pairs,1), numel(rest_files));
%         for n_pair = 1:size(pairs,1)
%             paris_all(n_pair, :) =  squeeze(corr_all_cat(pairs(n_pair,1),pairs(n_pair,2),rest_files));
%         end
%         [~,p] = ttest(paris_all(:,1),paris_all(:,end));
%         figure; hold on;
%         plot(labels_rest, paris_all, 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
%         errorbar(labels_rest, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
%         ylabel('pairwise correlations')
%         title(sprintf('rest pre-post pval=%.2f', p));
% 

        
        
    end
end

%% plot nobh ammn corr
corr_sig_ammn_nobh2 = cat(1,corr_sig_ammn_nobh{:});
corr_noise_ammn_nobh2 = cat(1,corr_noise_ammn_nobh{:});
rest_corr_all2  = cat(1,rest_corr_all{:});

[~,p] = ttest(corr_sig_ammn_nobh2(:,1),corr_sig_ammn_nobh2(:,end));
figure; hold on;
plot(sig_corr_labels, corr_sig_ammn_nobh2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(sig_corr_labels, mean(corr_sig_ammn_nobh2,1),std(corr_sig_ammn_nobh2,[],1)./sqrt(size(corr_sig_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('signal correlations ammn pre-post; pval=%.2f', p));


[~,p] = ttest(corr_noise_ammn_nobh2(:,1),corr_noise_ammn_nobh2(:,end));
figure; hold on;
plot(noise_corr_labels, corr_noise_ammn_nobh2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(noise_corr_labels, mean(corr_noise_ammn_nobh2,1),std(corr_noise_ammn_nobh2,[],1)./sqrt(size(corr_noise_ammn_nobh2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('noise correlations ammn pre-post; pval=%.2f', p));

labels_cell_rest = {'pre', 'post'};
labels_rest = reordercats(categorical(labels_cell_rest),labels_cell_rest);


[~,p] = ttest(rest_corr_all2(:,1),rest_corr_all2(:,end));
figure; hold on;
plot(labels_rest, rest_corr_all2', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(rest_corr_all2,1),std(rest_corr_all2,[],1)./sqrt(size(rest_corr_all2,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('rest pre-post pval=%.2f', p));

%%
% 
% 
% pairs = nchoosek(1:num_cells,2);
% 
% 
% labels_cell_ammn = {'pre', 'stim', 'post'};
% labels_ammn = reordercats(categorical(labels_cell_ammn),labels_cell_ammn);
% 
% SI_all_cat = cat(3,SI_all{:});
% paris_all = zeros(size(pairs,1), numel(ammn_files));
% for n_pair = 1:size(pairs,1)
%     paris_all(n_pair,:) =  squeeze(SI_all_cat(pairs(n_pair,1),pairs(n_pair,2),ammn_files));
% end
% [~,p] = ttest(paris_all(:,1),paris_all(:,end));
% figure; hold on;
% plot(labels_ammn, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
% errorbar(labels_ammn, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
% ylabel('cosine similarity')
% title(sprintf('ammn pre-post pval=%.2f', p));
% 
% 
% corr_all_cat = cat(3,corr_all{:});
% paris_all = zeros(size(pairs,1), numel(ammn_files));
% for n_pair = 1:size(pairs,1)
%     paris_all(n_pair, :) =  squeeze(corr_all_cat(pairs(n_pair,1),pairs(n_pair,2),ammn_files));
% end
% [~,p] = ttest(paris_all(:,1),paris_all(:,end));
% figure; hold on;
% plot(labels_ammn, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
% errorbar(labels_ammn, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
% ylabel('pairwise correlations')
% title(sprintf('ammn pre-post pval=%.2f', p));
% 

%%

fl_cell_traces_prsd2 = fl_cell_traces_prsd(ammn_stim_nobh);
fl_proc_data2 = fl_proc_data(ammn_stim_nobh);
fl_AC_data2 = fl_AC_data(ammn_stim_nobh);
fl_cell_type_stim2 = fl_cell_type_stim(ammn_stim_nobh);



rest_files = find(strcmpi([AC_data.paradigm], 'spont'));

pairs = nchoosek(1:num_cells,2);

labels_cell_rest = {'pre', 'post'};
labels_rest = reordercats(categorical(labels_cell_rest),labels_cell_rest);

SI_all_cat = cat(3,SI_all{:});
paris_all = zeros(size(pairs,1), numel(rest_files));
for n_pair = 1:size(pairs,1)
    paris_all(n_pair, :) =  squeeze(SI_all_cat(pairs(n_pair,1),pairs(n_pair,2),rest_files));
end
[~,p] = ttest(paris_all(:,1),paris_all(:,end));
figure; hold on;
plot(labels_rest, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('rest pre-post pval=%.2f', p));

corr_all_cat = cat(3,corr_all{:});
paris_all = zeros(size(pairs,1), numel(rest_files));
for n_pair = 1:size(pairs,1)
    paris_all(n_pair, :) =  squeeze(corr_all_cat(pairs(n_pair,1),pairs(n_pair,2),rest_files));
end
[~,p] = ttest(paris_all(:,1),paris_all(:,end));
figure; hold on;
plot(labels_rest, paris_all, 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('pairwise correlations')
title(sprintf('rest pre-post pval=%.2f', p));




%%

data = load_data.data;
params = load_data.params;
%%
AC_data = data.AC_data;
proc_data_all = data.proc_data;

ammn_files = find(strcmpi([AC_data.paradigm], 'ammn')+strcmpi([AC_data.paradigm], 'ammn_stim'));
rest_files = find(strcmpi([AC_data.paradigm], 'spont'));

labels_cell_ammn = {'pre', 'stim', 'post'};
labels_ammn = reordercats(categorical(labels_cell_ammn),labels_cell_ammn);

num_dsets = numel(AC_data.paradigm);


if sum(strcmpi([AC_data.paradigm], 'ammn_stim'))
    stat_idx = find(strcmpi([AC_data.paradigm], 'ammn_stim'));
elseif sum(strcmpi([AC_data.paradigm], 'spont_stim'))
    stat_idx = find(strcmpi([AC_data.paradigm], 'spont_stim'));
end
trace_type1 = cell(num_dsets,1);
for n_dset = 1:num_dsets
    if n_dset < stat_idx
        trace_type1{n_dset} = 'Pre';
    elseif n_dset > stat_idx
        trace_type1{n_dset} = 'Post';
    else
        trace_type1{n_dset} = 'Stim';
    end
end

trace_type = [trace_type1,[AC_data.paradigm]];

%%
% figure; hold on
% plot(proc_data_all{5, 1}.data.volt_data_all_aligned(:,1))
% plot(proc_data_all{5, 1}.data.volt_data_all_aligned(:,7))
% plot(proc_data_all{5, 1}.data.volt_data_all_aligned(:,8))
% 

traces_all = cat(1,data.cell_ca{:});
traces_alld = f_smooth_dfdt3(traces_all,1, 2, 1);


cell_run_id = cat(1,data.cell_run_id{:});
% 
% figure; imagesc(traces_alld(cell_run_id==1,:))
% figure; plot(mean(traces_alld(cell_run_id==1,:)))
%% Convert into cell with file in each
num_vol_all = zeros(num_dsets,1);
traces_alldspl = cell(num_dsets,1);
traces_allspl = cell(num_dsets,1);
cur_start = 1;
for n_dset = 1:num_dsets
    num_vol_all(n_dset) = proc_data_all{n_dset}.data.frame_data.num_volumes_linear;
    cur_end = cur_start +  num_vol_all(n_dset) - 1;
    traces_alldspl{n_dset} = traces_alld(:,cur_start:cur_end);
    traces_allspl{n_dset} = traces_all(:,cur_start:cur_end);
    cur_start = cur_end + 1;
end


%%

MMN_orientations = proc_data_all{ammn_files(1)}.data.MMN_orientations;
vol_period = proc_data_all{ammn_files(1)}.data.frame_data.volume_period;
plot_t = ((1:sum(base_onset_win))-base_onset_win(1))*vol_period/1000;

%%
cell_type_idx = cell_run_id==1;
num_cells = sum(cell_type_idx);

params.paradigm_to_check = {'ammn', 'ammn_stim'};
params.cell_type_idx = cell_type_idx;
params.corr_trials_check = {170};
params.trials_to_use_for_norm = [1:10 170 270];
%                     170, 270,...
%                      4, 7,... %,...
%                      201:206};
                     %1, 2, 3, 5, 6, 8, 9,...
                     %1 2 3 4 5 6 8 9 10,...
                     %201:206};%,...
                     %101:106};
params.base_onset_win = [0 15];

params.corr_method = 'cosine';
params.plot_sing_tr = 1;
params.plot_means = 1;
params.plot_corr = 1;
params.plot_rast = 0;
params.plot_sig_corr = 1;
params.norm_to_all_tr = 1;

params.set_mean_max = 0.7;
params.set_corr_ylim = [0 .5];
params.set_sig_corr_ylim = [0 .5];

f_s5_sig_noise_corr_within(traces_alldspl, proc_data_all, AC_data, params);

%%
                     

norm_to_all_tr = 1;

set_mean_max = 0.7;
set_corr_ylim = [0 .5];
set_sig_corr_ylim = [0 .5];

%set_corr_ylim = [-.1 .5];
%set_sig_corr_ylim = [-.1 .5];


norm_sig_corr_mat_all = cell(3,1);
norm_noise_corr_mat_all = cell(3,1);

for n_file_idx = 1:numel(ammn_files)
    n_file = ammn_files(n_file_idx);
    
    traces_all_d2 = traces_alldspl{n_file}(cell_type_idx,:);
    
    trial_types = proc_data_all{n_file}.data.trial_types;
    %stim_frame_index = proc_data_all{n_file}.data.stim_frame_index{1};
    
    stim_frame_index = proc_data_all{n_file}.data.stim_times_frame{1};
    
    stim_frame_index1 = stim_frame_index(logical(sum(trial_types == trials_to_use_for_norm,2)));
    trace_sort_temp = f_get_stim_trig_resp(traces_all_d2, stim_frame_index1, base_onset_win);
    
    % noise corr
    vec = reshape(trace_sort_temp, num_cells, [])';
    if strcmpi(corr_method, 'corr')
        vec = (vec - mean(vec,1));
    end
    norm_noise_corr_mat_all{n_file_idx} = (vecnorm(vec).*vecnorm(vec)');
    
    sig_all = cell(numel(trials_to_use_for_norm),1);
    for n_tt = 1:numel(trials_to_use_for_norm)
        stim_frame_index1 = stim_frame_index(logical(sum(trial_types == trials_to_use_for_norm(n_tt),2)));
        temp_trace = f_get_stim_trig_resp(traces_all_d2, stim_frame_index1, base_onset_win);
        sig_all{n_tt} = mean(temp_trace,3);
    end
    vec = cat(2,sig_all{:})';
    if strcmpi(corr_method, 'corr')
        vec = (vec - mean(vec,1));
    end
    norm_sig_corr_mat_all{n_file_idx} = (vecnorm(vec).*vecnorm(vec)');
end


for n_tt = 1:numel(corr_trials_check)
    tt1 = corr_trials_check{n_tt};
    
    trace_sort_all = cell(numel(ammn_files), 1);
    trial_types_all = cell(numel(ammn_files), 1);
    for n_file_idx = 1:numel(ammn_files)
        n_file = ammn_files(n_file_idx);
        traces_all_d2 = traces_alldspl{n_file}(cell_type_idx,:);
        trial_types = proc_data_all{n_file}.data.trial_types;
        %stim_frame_index = proc_data_all{n_file}.data.stim_frame_index{1};
        stim_frame_index = proc_data_all{n_file}.data.stim_times_frame{1,1};
        stim_frame_index1 = stim_frame_index(logical(sum(trial_types == tt1,2)));
        num_trials = numel(stim_frame_index1);

        trace_sort_all{n_file_idx} = f_get_stim_trig_resp(traces_all_d2, stim_frame_index1, base_onset_win);
    end
    
    % compute normalization mat
    trials_to_use_for_norm = [1:10 170 270];

    
    if plot_sing_tr
        % some single cell resp
        for n_file_idx = 1:numel(ammn_files)
            n_file = ammn_files(n_file_idx);
            trace_sort = trace_sort_all{n_file_idx};
            figure; 
            for n_cell = 1:num_cells
                subplot(1, num_cells, n_cell); hold on; axis tight;
                plot(plot_t, squeeze(trace_sort(n_cell,:,:)), 'color', [.6 .6 .6])
                plot(plot_t, squeeze(mean(trace_sort(n_cell,:,:),3)), 'm', 'linewidth', 2)
                title(sprintf('cell %d', n_cell))
                ylim([0 1])
                if n_cell == 1
                    xlabel('Time (sec)')
                end
            end
            sgtitle(sprintf('%s %s trial %s', trace_type{n_file,2}, trace_type{n_file,1}, num2str(tt1)))
        end
    end
    
    if plot_means
        if set_mean_max
            y_lim2 = set_mean_max;
        else
            y_lim2 = 0.01;
            for n_file_idx = 1:numel(ammn_files)
                y_lim2 = max([y_lim2 max(max(mean(trace_sort_all{n_file_idx},3)))]);
            end
        end
        stim_colors = {[77 206 124]/255, [245 99 119]/255, [21 69 225]/255};

        figure;
        for n_cell = 1:num_cells
            subplot(1, num_cells, n_cell); hold on; axis tight;
            for n_file_idx = 1:numel(ammn_files)
                trace_sort = trace_sort_all{n_file_idx};
                plot(plot_t, squeeze(mean(trace_sort(n_cell,:,:),3)), 'linewidth', 2, 'color', stim_colors{n_file_idx})
            end
            ylim([0 y_lim2])
            if n_cell == 1
                xlabel('Time (sec)')
            end
        end
        legend('pre', 'stim', 'post')
        sgtitle(sprintf('Mean response to trial %s', num2str(tt1)))
    end
    
    
    if strcmpi(corr_method, 'corr')
        title_tag = 'pairwise correlations';
    elseif strcmpi(corr_method, 'cosine')
        title_tag = 'cosine similarity';
    end
    
    if plot_corr
        % pairwise corr noise
        corr_tt = cell(numel(ammn_files), 1);
        
        if norm_to_all_tr
            for n_file_idx = 1:numel(ammn_files)
                vec = reshape(trace_sort_all{n_file_idx}, num_cells, [])';
                if strcmpi(corr_method, 'corr')
                    vec = (vec - mean(vec,1));
                end
                corr_tt{n_file_idx} = (vec'*vec)./norm_noise_corr_mat_all{n_file_idx};
            end
        else
            if strcmpi(corr_method, 'corr')
                for n_file_idx = 1:numel(ammn_files)
                    corr_tt{n_file_idx} = corr(reshape(trace_sort_all{n_file_idx}, num_cells, [])');
                end
            elseif strcmpi(corr_method, 'cosine')
                for n_file_idx = 1:numel(ammn_files)
                    corr_tt{n_file_idx} = 1-f_pdist_YS(reshape(trace_sort_all{n_file_idx}, num_cells, []), 'cosine');
                end
            end
        end
        
        
        pairs = nchoosek(1:num_cells,2);

        corr_tt_cat = cat(3,corr_tt{:});
        paris_all = zeros(size(pairs,1), numel(ammn_files));
        for n_pair = 1:size(pairs,1)
            paris_all(n_pair,:) =  squeeze(corr_tt_cat(pairs(n_pair,1),pairs(n_pair,2),:));
        end
        [~,p] = ttest(paris_all(:,1),paris_all(:,end));
        figure; hold on;
        plot(labels_ammn, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
        errorbar(labels_ammn, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
        ylabel(title_tag)
        title(sprintf('noise correlations ammn pre-post; trial=%s; pval=%.2f', num2str(tt1), p));
        if ~isempty(set_corr_ylim)
            ylim(set_corr_ylim)
        end
    end
    
    if plot_rast
        for n_file_idx = 1:numel(ammn_files)
            n_file = ammn_files(n_file_idx);
            temp_trace = reshape(trace_sort_all{n_file_idx}, num_cells, []);
            figure; imagesc(temp_trace);
            caxis([0 .7]);
            title(sprintf('%s %s trial %s', trace_type{n_file,2}, trace_type{n_file,1}, num2str(tt1)));
        end

        figure; hold on
        for n_file_idx = 1:numel(ammn_files)
            n_file = ammn_files(n_file_idx);
            temp_trace = reshape(trace_sort_all{n_file_idx}, num_cells, []);

            temp_mat = temp_trace*temp_trace';
            %mean([temp_mat(2,1), temp_mat(3,1), temp_mat(3,2)]);
            plot(temp_trace(1,:).*temp_trace(2,:).*temp_trace(3,:))
        end
    end
    
    if plot_sig_corr
        % signal correlations
        
        corr_sig_tt = cell(numel(ammn_files), 1);
        
        if norm_to_all_tr
            for n_file_idx = 1:numel(ammn_files)
                vec = mean(trace_sort_all{n_file_idx},3)';
                if strcmpi(corr_method, 'corr')
                    vec = (vec - mean(vec,1));
                end
                corr_sig_tt{n_file_idx} = (vec'*vec)./norm_sig_corr_mat_all{n_file_idx};
            end
        else
            if strcmpi(corr_method, 'corr')
                for n_file_idx = 1:numel(ammn_files)
                    corr_sig_tt{n_file_idx} = corr(mean(trace_sort_all{n_file_idx},3)');
                end
            elseif strcmpi(corr_method, 'cosine')
                for n_file_idx = 1:numel(ammn_files)
                    corr_sig_tt{n_file_idx} = 1-f_pdist_YS(mean(trace_sort_all{n_file_idx},3), 'cosine');
                end
            end
        end
        
        pairs = nchoosek(1:num_cells,2);

        corr_tt_cat = cat(3,corr_sig_tt{:});
        paris_all = zeros(size(pairs,1), numel(ammn_files));
        for n_pair = 1:size(pairs,1)
            paris_all(n_pair,:) =  squeeze(corr_tt_cat(pairs(n_pair,1),pairs(n_pair,2),:));
        end
        [~,p] = ttest(paris_all(:,1),paris_all(:,end));
        figure; hold on;
        plot(labels_ammn, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
        errorbar(labels_ammn, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
        ylabel(title_tag)
        title(sprintf('signal correlations ammn pre-post; trial=%s; pval=%.2f', num2str(tt1), p));
        if ~isempty(set_sig_corr_ylim)
            ylim(set_sig_corr_ylim)
        end
        

    end
end

%%

% 
% % signal correlations 2 ... normalize by whole trace correlation
% corr_sig_tt = cell(numel(ammn_files), 1);
% if strcmpi(corr_method, 'corr')
%     for n_file_idx = 1:numel(ammn_files)
%         corr_sig_tt{n_file_idx} = corr(mean(trace_sort_all{n_file_idx},3)');
%         title_tag = 'pairwise correlations';
%     end
% elseif strcmpi(corr_method, 'cosine')
%     for n_file_idx = 1:numel(ammn_files)
%         corr_sig_tt{n_file_idx} = 1-f_pdist_YS(mean(trace_sort_all{n_file_idx},3), 'cosine');
%         title_tag = 'cosine similarity';
%     end
% end
% pairs = nchoosek(1:num_cells,2);
% 
% corr_tt_cat = cat(3,corr_sig_tt{:});
% paris_all = zeros(size(pairs,1), numel(ammn_files));
% for n_pair = 1:size(pairs,1)
%     paris_all(n_pair,:) =  squeeze(corr_tt_cat(pairs(n_pair,1),pairs(n_pair,2),:));
% end
% [~,p] = ttest(paris_all(:,1),paris_all(:,end));
% figure; hold on;
% plot(labels_ammn, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
% errorbar(labels_ammn, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
% ylabel(title_tag)
% title(sprintf('ammn pre-post signal correlations; trial=%s; pval=%.2f', num2str(tt1), p));
% if ~isempty(set_sig_corr_ylim)
%     ylim(set_sig_corr_ylim)
% end




%%

if num_dsets == 5
    plot_locs = [4 1 2 6 3];
    spm = 2;
    spn = 3;
elseif num_dsets == 6
    plot_locs = [1 2 6 3 8 4];
    spm = 2;
    spn = 4;
end

SI_all = cell(num_dsets, 1);

for n_file = 1:num_dsets
    SI_all{n_file} = 1-f_pdist_YS(traces_alldspl{n_file}(cell_type_idx,:), 'cosine');
end

SI_mean = mean(cat(3,SI_all{:}),3);

SI_all_n = SI_all;
for n_file = 1:num_dsets
    SI_all_n{n_file} = SI_all{n_file} - SI_mean;
end

clim1 = [min(min(min(cat(3,SI_all{:})))) max(max(max(cat(3,SI_all{:}))))];
clim2 = [min(min(min(cat(3,SI_all_n{:})))) max(max(max(cat(3,SI_all_n{:}))))];

clim3 = [.3 .6];


figure;
for n_file = 1:num_dsets
    sp = subplot(spm,spn, plot_locs(n_file));
    imagesc(SI_all{n_file})
    caxis(clim3);
    %sp.XTick = [];
    %sp.YTick = [];
    title(sprintf('%s %s', trace_type{n_file,2}, trace_type{n_file,1}), 'interpreter', 'none')
    if plot_locs(n_file) == 1
        ylabel('Cells')
    end
    if plot_locs(n_file) == 4
        ylabel('Cells')
        xlabel('Cells')
    end
end
sgtitle('cell-cell cosine similarity')



%%
corr_all = cell(num_dsets, 1);

for n_file = 1:num_dsets
    corr_all{n_file} = corr(traces_alldspl{n_file}(cell_type_idx,:)');
end

figure;
for n_file = 1:num_dsets
    sp = subplot(spm,spn, plot_locs(n_file));
    imagesc(corr_all{n_file})
    caxis([0 .5]);
    %sp.XTick = [];
    %sp.YTick = [];
    title(sprintf('%s %s', trace_type{n_file,2}, trace_type{n_file,1}))
    if plot_locs(n_file) == 1
        ylabel('Cells')
    end
    if plot_locs(n_file) == 4
        ylabel('Cells')
        xlabel('Cells')
    end
end
sgtitle('pairwise corr')

%%

pairs = nchoosek(1:num_cells,2);


labels_cell_ammn = {'pre', 'stim', 'post'};
labels_ammn = reordercats(categorical(labels_cell_ammn),labels_cell_ammn);

SI_all_cat = cat(3,SI_all{:});
paris_all = zeros(size(pairs,1), numel(ammn_files));
for n_pair = 1:size(pairs,1)
    paris_all(n_pair,:) =  squeeze(SI_all_cat(pairs(n_pair,1),pairs(n_pair,2),ammn_files));
end
[~,p] = ttest(paris_all(:,1),paris_all(:,end));
figure; hold on;
plot(labels_ammn, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_ammn, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('ammn pre-post pval=%.2f', p));


corr_all_cat = cat(3,corr_all{:});
paris_all = zeros(size(pairs,1), numel(ammn_files));
for n_pair = 1:size(pairs,1)
    paris_all(n_pair, :) =  squeeze(corr_all_cat(pairs(n_pair,1),pairs(n_pair,2),ammn_files));
end
[~,p] = ttest(paris_all(:,1),paris_all(:,end));
figure; hold on;
plot(labels_ammn, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_ammn, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('pairwise correlations')
title(sprintf('ammn pre-post pval=%.2f', p));



labels_cell_rest = {'pre', 'post'};
labels_rest = reordercats(categorical(labels_cell_rest),labels_cell_rest);

SI_all_cat = cat(3,SI_all{:});
paris_all = zeros(size(pairs,1), numel(rest_files));
for n_pair = 1:size(pairs,1)
    paris_all(n_pair, :) =  squeeze(SI_all_cat(pairs(n_pair,1),pairs(n_pair,2),rest_files));
end
[~,p] = ttest(paris_all(:,1),paris_all(:,end));
figure; hold on;
plot(labels_rest, paris_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('cosine similarity')
title(sprintf('rest pre-post pval=%.2f', p));

corr_all_cat = cat(3,corr_all{:});
paris_all = zeros(size(pairs,1), numel(rest_files));
for n_pair = 1:size(pairs,1)
    paris_all(n_pair, :) =  squeeze(corr_all_cat(pairs(n_pair,1),pairs(n_pair,2),rest_files));
end
[~,p] = ttest(paris_all(:,1),paris_all(:,end));
figure; hold on;
plot(labels_rest, paris_all, 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
errorbar(labels_rest, mean(paris_all,1),std(paris_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
ylabel('pairwise correlations')
title(sprintf('rest pre-post pval=%.2f', p));


%% 

figure;
for n_file = 1:num_files
    sp = subplot(spm,3, plot_locs(n_file));
    imagesc(SI_all_n{n_file})
    caxis(clim2);
    %sp.XTick = [];
    %sp.YTick = [];
    title(sprintf('%s %s', trace_type{n_file,2}, trace_type{n_file,1}))
    if plot_locs(n_file) == 1
        ylabel('Cells')
    end
    if plot_locs(n_file) == 4
        ylabel('Cells')
        xlabel('Cells')
    end
end
sgtitle('cell-cell cosine similarity')




%%
if 1
    % here need to do stim trig ave
    if sum(strcmpi([AC_data.paradigm], 'behavior'))
        bh_proc = proc_data_all{strcmpi([AC_data.paradigm], 'behavior')};
        bh_freq_ops = bh_proc.data.stim_params.ops;
        stim_tone_bh = bh_freq_ops.dev_tone_list;
    end
    
    if sum(strcmpi([AC_data.paradigm], 'ammn_stim'))
        ammn_stim_proc = proc_data_all{strcmpi([AC_data.paradigm], 'ammn_stim')};
        ammn_freq_ops = ammn_stim_proc.data.stim_params.ops;
        stim_run = ammn_freq_ops.stim_trials_volt{1};
        mmn_pat = ammn_freq_ops.MMN_patterns(ammn_freq_ops.paradigm_MMN_pattern(stim_run),:);
        stim_tone_dd = (stim_run - 1)*100 + 70;
        stim_tone_cont = mmn_pat(2);
    end
    
    num_frames_cumsum = cumsum(num_vol_all);
    win1 = [5 60];
    
    trials_check = [5, 170];
    
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
        
        if sum(strcmpi(AC_data.paradigm{n_pr},{'behavior', 'ammn_stim', 'spont_stim'}))
            if sum(strcmpi(AC_data.paradigm{n_pr}, {'behavior'}))
                stim_chan = 'stim type';
            elseif sum(strcmpi(AC_data.paradigm{n_pr}, {'ammn_stim', 'spont_stim'}))
                stim_chan = 'pockel';
            end
            idx_chan = strcmpi(proc_data_all{n_pr}.ops.chan_labels, stim_chan);
            stim_times_frame = proc_data_all{n_pr}.data.stim_times_frame{idx_chan,1};
            
            stim_times_frame(stim_times_frame < win1(1)) = [];
            stim_times_frame((stim_times_frame + win1(2)) > num_vol_all(n_pr)) = [];
            
            trace1n = traces_allspl{n_pr}(cell_type_idx,:);
            trace1nd = traces_alldspl{n_pr}(cell_type_idx,:);
            %stim_trace2 = trace1n((num_frames_cumsum(n_pr)-num_vol_all(n_pr)+1):num_frames_cumsum(n_pr));
            sorted_trials = f_get_stim_trig_resp(trace1n, stim_times_frame, win1);
            %tim_trace2d = trace1nd((num_frames_cumsum(n_pr)-num_vol_all(n_pr)+1):num_frames_cumsum(n_pr));
            sorted_trialsd = f_get_stim_trig_resp(trace1nd, stim_times_frame, win1);
            
            for n_cell = 1:size(trace1n,1)
                cell_sorted_trials = squeeze(sorted_trials(n_cell,:,:));
                cell_sorted_trialsd = squeeze(sorted_trialsd(n_cell,:,:));
                figure; 
                subplot(3,1,1);
                imagesc(cell_sorted_trials'); 
                title('trials all');
                subplot(3,1,2); hold on;
                plot(cell_sorted_trials, 'color', [.6 .6 .6])
                plot(mean(cell_sorted_trials,2), 'Linewidth', 2, 'color', 'k');
                axis tight;
                title('trig ave');
                subplot(3,1,3); hold on;
                plot(cell_sorted_trialsd, 'color', [.6 .6 .6])
                plot(mean(cell_sorted_trialsd,2), 'Linewidth', 2, 'color', 'k');
                axis tight;
                title('trig ave deconvolved');
                sgtitle(sprintf('%s; cell %d, %s',stim_chan, n_cell, AC_data.paradigm{n_pr}));
            end
        end
        
        if sum(strcmpi(AC_data.paradigm{n_pr}, {'ammn', 'ammn_stim'}))
            stim_chan = 'stim type';
            idx_chan = strcmpi(proc_data_all{n_pr}.ops.chan_labels, stim_chan);
            
            trial_types = proc_data_all{n_pr}.data.trial_types;
            
            for n_tr = 1:numel(trials_check)
                trial_idx = (trial_types == trials_check(n_tr));
                stim_times_frame = proc_data_all{n_pr}.data.stim_times_frame{idx_chan,1};
                stim_times_frame = stim_times_frame(trial_idx);
                stim_times_frame(stim_times_frame < win1(1)) = [];
                stim_times_frame((stim_times_frame + win1(2)) > num_vol_all(n_pr)) = [];

                trace1n = traces_allspl{n_pr}(cell_type_idx,:);
                trace1nd = traces_alldspl{n_pr}(cell_type_idx,:);
                %stim_trace2 = trace1n((num_frames_cumsum(n_pr)-num_vol_all(n_pr)+1):num_frames_cumsum(n_pr));
                sorted_trials = f_get_stim_trig_resp(trace1n, stim_times_frame, win1);
                %tim_trace2d = trace1nd((num_frames_cumsum(n_pr)-num_vol_all(n_pr)+1):num_frames_cumsum(n_pr));
                sorted_trialsd = f_get_stim_trig_resp(trace1nd, stim_times_frame, win1);

                for n_cell = 1:size(trace1n,1)
                    cell_sorted_trials = squeeze(sorted_trials(n_cell,:,:));
                    cell_sorted_trialsd = squeeze(sorted_trialsd(n_cell,:,:));
                    figure; 
                    subplot(3,1,1);
                    imagesc(cell_sorted_trials'); 
                    title('trials all');
                    subplot(3,1,2); hold on;
                    plot(cell_sorted_trials, 'color', [.6 .6 .6])
                    plot(mean(cell_sorted_trials,2), 'Linewidth', 2, 'color', 'k');
                    axis tight;
                    title('trig ave');
                    subplot(3,1,3); hold on;
                    plot(cell_sorted_trialsd, 'color', [.6 .6 .6])
                    plot(mean(cell_sorted_trialsd,2), 'Linewidth', 2, 'color', 'k');
                    axis tight;
                    title('trig ave deconvolved');
                    sgtitle(sprintf('%s; cell %d; trial %d; %s',stim_chan, n_cell, trials_check(n_tr), AC_data.paradigm{n_pr}));
                end
            end
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