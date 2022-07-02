function data_out = f_s5_sig_noise_corr_within(traces_alldspl, proc_data_all, AC_data, params)

corr_method = 'cosine';

plot_sing_tr = 1;
plot_means = 1;
plot_corr = 1;
plot_rast = 0;
plot_sig_corr = 1;

paradigm_to_check = params.paradigm_to_check;

ammn_files = [];
for n_par = 1:numel(paradigm_to_check)
    x = find(strcmpi([AC_data.paradigm], paradigm_to_check{n_par}));
    ammn_files = [ammn_files; x];
end
ammn_files = unique(ammn_files);

base_onset_win = params.base_onset_win;

MMN_orientations = proc_data_all{ammn_files(1)}.data.MMN_orientations;
vol_period = proc_data_all{ammn_files(1)}.data.frame_data.volume_period;
plot_t = ((1:sum(base_onset_win))-base_onset_win(1))*vol_period/1000;


cell_type_idx = params.cell_type_idx;
corr_trials_check = params.corr_trials_check;
%corr_trials_check = {170};
%                     170, 270,...
%                      4, 7,... %,...
%                      201:206};
                     %1, 2, 3, 5, 6, 8, 9,...
                     %1 2 3 4 5 6 8 9 10,...
                     %201:206};%,...
                     %101:106};
                     
trials_to_use_for_norm = params.trials_to_use_for_norm;

norm_to_all_tr = params.norm_to_all_tr;

set_mean_max = params.set_mean_max;
set_corr_ylim = params.set_corr_ylim;
set_sig_corr_ylim = params.set_sig_corr_ylim;

num_dsets = numel(AC_data.paradigm);
if sum(strcmpi([AC_data.paradigm], 'ammn_stim'))
    stim_idx = find(strcmpi([AC_data.paradigm], 'ammn_stim'));
elseif sum(strcmpi([AC_data.paradigm], 'spont_stim'))
    stim_idx = find(strcmpi([AC_data.paradigm], 'spont_stim'));
end
trace_type1 = cell(num_dsets,1);
for n_dset = 1:num_dsets
    if n_dset < stim_idx
        trace_type1{n_dset} = 'Pre';
    elseif n_dset > stim_idx
        trace_type1{n_dset} = 'Post';
    else
        trace_type1{n_dset} = 'Stim';
    end
end
trace_type = [trace_type1,[AC_data.paradigm]];

labels_cell_ammn = {'pre', 'stim', 'post'};
labels_ammn = reordercats(categorical(labels_cell_ammn),labels_cell_ammn);


norm_sig_corr_mat_all = cell(3,1);
norm_noise_corr_mat_all = cell(3,1);

num_cells = sum(cell_type_idx);
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

    
    if params.plot_sing_tr
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
    
    if params.plot_means
        if set_mean_max
            y_lim2 = set_mean_max;
        else
            y_lim2 = 0.01;
            for n_file_idx = 1:numel(ammn_files)
                y_lim2 = max([y_lim2 max(max(mean(trace_sort_all{n_file_idx},3)))]);
            end
        end
        stim_colors = {[77 206 124]/255, [245 99 119]/255, [21 69 225]/255, [1 1 1]};

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
    
    if params.plot_corr
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
        pairs_all = zeros(size(pairs,1), numel(ammn_files));
        for n_pair = 1:size(pairs,1)
            pairs_all(n_pair,:) =  squeeze(corr_tt_cat(pairs(n_pair,1),pairs(n_pair,2),:));
        end
        [~,p] = ttest(pairs_all(:,1),pairs_all(:,end));
        figure; hold on;
        plot(labels_ammn, pairs_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
        errorbar(labels_ammn, mean(pairs_all,1),std(pairs_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
        ylabel(title_tag)
        title(sprintf('noise correlations ammn pre-post; trial=%s; pval=%.2f', num2str(tt1), p));
        if ~isempty(set_corr_ylim)
            ylim(set_corr_ylim)
        end
        
        data_out.noise_corr.labels_ammn = labels_ammn;
        data_out.noise_corr.pairs_all = pairs_all;
    end
    
    if params.plot_rast
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
    
    if params.plot_sig_corr
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
        pairs_all = zeros(size(pairs,1), numel(ammn_files));
        for n_pair = 1:size(pairs,1)
            pairs_all(n_pair,:) =  squeeze(corr_tt_cat(pairs(n_pair,1),pairs(n_pair,2),:));
        end
        [~,p] = ttest(pairs_all(:,1),pairs_all(:,end));
        figure; hold on;
        plot(labels_ammn, pairs_all', 'o-', 'Linewidth', 1, 'color', [.6 .6 .6])
        errorbar(labels_ammn, mean(pairs_all,1),std(pairs_all,[],1)./sqrt(size(pairs,1)-1), 'o-', 'Linewidth', 2, 'color', [.2 .2 .2])
        ylabel(title_tag)
        title(sprintf('signal correlations ammn pre-post; trial=%s; pval=%.2f', num2str(tt1), p));
        if ~isempty(set_sig_corr_ylim)
            ylim(set_sig_corr_ylim)
        end
        
        data_out.sig_corr.labels_ammn = labels_ammn;
        data_out.sig_corr.pairs_all = pairs_all;

    end
end




end