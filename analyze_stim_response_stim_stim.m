

clear;
close all;

fpath_preproc = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data_dream\preprocessing\';
fpath_cells = 'D:\data\AC\M4460_1_2_22b_dream\cell_traces\';

base_onset_win = [0 9];

%%

raw_traces_list = {'cell1_ammn1.csv', 'cell2_ammn1.csv', 'cell3_ammn1.csv', 'cell4_ammn1.csv', 'cell5_ammn1.csv', 'cell6_ammn1.csv';...
                   'cell1_rest2.csv', 'cell2_rest2.csv', 'cell3_rest2.csv', 'cell4_rest2.csv', 'cell5_rest2.csv', 'cell6_rest2.csv';...
                   'cell1_ammn_stim3.csv', 'cell2_ammn_stim3.csv', 'cell3_ammn_stim3.csv', 'cell4_ammn_stim3.csv', 'cell5_ammn_stim3.csv', 'cell6_ammn_stim3.csv';...
                   'cell1_rest4.csv', 'cell2_rest4.csv', 'cell3_rest4.csv', 'cell4_rest4.csv', 'cell5_rest4.csv', 'cell6_rest4.csv';...
                   'cell1_ammn5.csv', 'cell2_ammn5.csv', 'cell3_ammn5.csv', 'cell4_ammn5.csv', 'cell5_ammn5.csv', 'cell6_ammn5.csv'};
               
trace_type = {'pre', 'ammn';...
              'pre', 'rest';...
              'stim', 'ammn';...
              'post', 'rest';...
              'post', 'ammn'};
cuts_stim_fname = {'M4460_im1_AC_ammn1_1_2_22b_h5cutsdata.mat';...
                   'M4460_im2_AC_rest2_1_2_22b_h5cutsdata.mat';...
                   'M4460_im3_AC_ammn_stim3_1_2_22b_h5cutsdata.mat';...
                   'M4460_im4_AC_rest4_1_2_22b_h5cutsdata.mat';...
                   'M4460_im5_AC_ammn5_1_2_22b_h5cutsdata.mat'};
               
preprocessed_data = {'M4460_im1_AC_ammn1_1_2_22b_processed_data.mat';
                     '';...
                     'M4460_im3_AC_ammn_stim3_1_2_22b_processed_data.mat';...
                     '';...
                     'M4460_im5_AC_ammn5_1_2_22b_processed_data.mat'};

%% load traces
[num_files, num_cells] = size(raw_traces_list);

traces_all = cell(num_files, 1);
traces_all_s = cell(num_files, 1);
traces_all_d = cell(num_files, 1);
cuts_all = cell(num_files, 1);
proc_data_all = cell(num_files, 1);

for n_file = 1:num_files
    temp_cuts = load([fpath_preproc '\' cuts_stim_fname{n_file}]);
    cuts_all{n_file} = temp_cuts.cuts_data{1};

    traces2 = cell(num_cells, 1);
    for n_cell = 1:num_cells
        trace_temp = csvread([fpath_cells '/' raw_traces_list{n_file, n_cell}],1,1);
        trace_full_temp = zeros(1, numel(cuts_all{n_file}.vid_cuts_trace));
        trace_full_temp(logical(cuts_all{n_file}.vid_cuts_trace)) = trace_temp;
        traces2{n_cell} = trace_full_temp;
    end
    traces_full = cat(1,traces2{:});
    traces_all{n_file} = traces_full;
    traces_all_d{n_file} = f_smooth_dfdt3(traces_full, 1, 2, 1);
    
    traces_all_s{n_file} = f_smooth_gauss2(traces_full, 4, 1);
    
    if ~isempty(preprocessed_data{n_file})
        proc_data_all{n_file} = load([fpath_preproc preprocessed_data{n_file}]);
    end
end

ammn_files = find(strcmpi(trace_type(:,2), 'ammn'));
rest_files = find(strcmpi(trace_type(:,2), 'rest'));


labels_cell_ammn = {'pre', 'stim', 'post'};
labels_ammn = reordercats(categorical(labels_cell_ammn),labels_cell_ammn);


%%

figure; hold on;
plot(traces_all_d{4}')

%%

MMN_orientations = proc_data_all{ammn_files(1)}.data.MMN_orientations;
vol_period = proc_data_all{ammn_files(1)}.data.frame_data.volume_period;
plot_t = ((1:sum(base_onset_win))-base_onset_win(1))*vol_period/1000;

%%

corr_method = 'cosine';

plot_sing_tr = 1;
plot_means = 1;
plot_corr = 1;
plot_rast = 0;
plot_sig_corr = 1;

corr_trials_check = {4, 270};
%                     170, 270,...
%                      4, 7,... %,...
%                      201:206};
                     %1, 2, 3, 5, 6, 8, 9,...
                     %1 2 3 4 5 6 8 9 10,...
                     %201:206};%,...
                     %101:106};
                     
trials_to_use_for_norm = [1:10 170 270];
norm_to_all_tr = 1;

set_mean_max = 0.7;
set_corr_ylim = [0 1];
set_sig_corr_ylim = [0 1];

%set_corr_ylim = [-.1 .5];
%set_sig_corr_ylim = [-.1 .5];

norm_sig_corr_mat_all = cell(3,1);
norm_noise_corr_mat_all = cell(3,1);
for n_file_idx = 1:numel(ammn_files)
    n_file = ammn_files(n_file_idx);
    trial_types = proc_data_all{n_file}.data.trial_types;
    stim_frame_index = proc_data_all{n_file}.data.stim_frame_index{1};
    stim_frame_index1 = stim_frame_index(logical(sum(trial_types == trials_to_use_for_norm,2)));
    trace_sort_temp = f_get_stim_trig_resp(traces_all_d{n_file}, stim_frame_index1, base_onset_win);
    
    % noise corr
    vec = reshape(trace_sort_temp, num_cells, [])';
    if strcmpi(corr_method, 'corr')
        vec = (vec - mean(vec,1));
    end
    norm_noise_corr_mat_all{n_file_idx} = (vecnorm(vec).*vecnorm(vec)');
    
    sig_all = cell(numel(trials_to_use_for_norm),1);
    for n_tt = 1:numel(trials_to_use_for_norm)
        stim_frame_index1 = stim_frame_index(logical(sum(trial_types == trials_to_use_for_norm(n_tt),2)));
        temp_trace = f_get_stim_trig_resp(traces_all_d{n_file}, stim_frame_index1, base_onset_win);
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
        trial_types = proc_data_all{n_file}.data.trial_types;
        stim_frame_index = proc_data_all{n_file}.data.stim_frame_index{1};
        stim_frame_index1 = stim_frame_index(logical(sum(trial_types == tt1,2)));
        num_trials = numel(stim_frame_index1);

        trace_sort_all{n_file_idx} = f_get_stim_trig_resp(traces_all_d{n_file}, stim_frame_index1, base_onset_win);
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


        % signal correlations 2 ... normalize by whole trace correlation
        corr_sig_tt = cell(numel(ammn_files), 1);
        if strcmpi(corr_method, 'corr')
            for n_file_idx = 1:numel(ammn_files)
                corr_sig_tt{n_file_idx} = corr(mean(trace_sort_all{n_file_idx},3)');
                title_tag = 'pairwise correlations';
            end
        elseif strcmpi(corr_method, 'cosine')
            for n_file_idx = 1:numel(ammn_files)
                corr_sig_tt{n_file_idx} = 1-f_pdist_YS(mean(trace_sort_all{n_file_idx},3), 'cosine');
                title_tag = 'cosine similarity';
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
        title(sprintf('ammn pre-post signal correlations; trial=%s; pval=%.2f', num2str(tt1), p));
        if ~isempty(set_sig_corr_ylim)
            ylim(set_sig_corr_ylim)
        end
        



%%

SI_all = cell(num_files, 1);

for n_file = 1:num_files
    SI_all{n_file} = 1-f_pdist_YS(traces_all_d{n_file}, 'cosine');
end

SI_mean = mean(cat(3,SI_all{:}),3);

SI_all_n = SI_all;
for n_file = 1:num_files
    SI_all_n{n_file} = SI_all{n_file} - SI_mean;
end

clim1 = [min(min(min(cat(3,SI_all{:})))) max(max(max(cat(3,SI_all{:}))))];
clim2 = [min(min(min(cat(3,SI_all_n{:})))) max(max(max(cat(3,SI_all_n{:}))))];

clim3 = [.3 .6];
plot_locs = [4 1 2 6 3];

figure;
for n_file = 1:num_files
    sp = subplot(2,3, plot_locs(n_file));
    imagesc(SI_all{n_file})
    caxis(clim3);
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
corr_all = cell(num_files, 1);

for n_file = 1:num_files
    corr_all{n_file} = corr(traces_all_d{n_file}');
end

figure;
for n_file = 1:num_files
    sp = subplot(2,3, plot_locs(n_file));
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
    sp = subplot(2,3, plot_locs(n_file));
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

