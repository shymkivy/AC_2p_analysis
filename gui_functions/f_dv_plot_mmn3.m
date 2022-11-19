function f_dv_plot_mmn3(app)

add_combined = 1;

n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = numel(data.experiment);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

num_t = sum(trial_window);

ctx_plot_list = [18, 19, 20; ...
                     28, 29, 30]';
                 
num_flip = size(ctx_plot_list,2);
                 
params = f_dv_gather_params(app);


reg_data = cell(4,1);
for n_reg = 1:4
    reg_data{n_reg}.resp_all = cell(num_dsets, num_flip);
    reg_data{n_reg}.cell_counts = zeros(num_dsets, num_flip);
end
for n_flip = 1:num_flip
    tn_all = ctx_plot_list(:,n_flip);
    ctx1 = app.ops.context_types_all(tn_all)';

    for n_dset = 1:num_dsets
        fprintf('Dset %d\n', n_dset);
        data1 =  data(n_dset,:);
        stats1 = data1.stats{n_pl};
        params.n_dset = find(data1.idx == app.data.idx);

        cdata = f_dv_compute_cdata(app.ddata, params);

        firing_rate = cdata.S_sm;
        trial_types = data1.trial_types{1};
        stim_times = data1.stim_frame_index{n_pl};
        mmn_freq = data1.MMN_freq{1};
        cell_is_resp = stats1.peak_resp_cells;

        if ~isempty(data1.registered_data{1})
            reg_labels = data1.registered_data{1}.reg_labels;
        else
            reg_labels = zeros(cdata.num_cells,1);
        end

        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
        [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

        if app.ConverttoZCheckBox.Value
            st_mean_mean = stats1.stat_trials_mean_mean;
            st_mean_sem = stats1.stat_trials_mean_sem;
        else
            st_mean_mean = zeros(cdata.num_cells,1);
            st_mean_sem = ones(cdata.num_cells,1);
        end
        
        trial_data_sort_wctx = (trial_data_sort_wctx - st_mean_mean)./st_mean_sem;
        
        resp_cell_idx = logical(sum(cell_is_resp(:,tn_all),2));
        for n_reg = 1:4
            reg_idx = reg_labels == n_reg;
            % get resp cells
            resp_cell_idx2 = logical(resp_cell_idx.*reg_idx);
            
            reg_data{n_reg}.cell_counts(n_dset, n_flip) = sum(resp_cell_idx2);

            reg_data{n_reg}.resp_all{n_dset, n_flip} = zeros(reg_data{n_reg}.cell_counts(n_dset, n_flip), num_t, size(ctx1,2));
            for n_ctx = 1:size(ctx1,2)
                ctx2 = ctx1(:,n_ctx)';
                temp_resp = trial_data_sort_wctx(resp_cell_idx2,:,logical(sum(trial_types_wctx == ctx2,2)));
                reg_data{n_reg}.resp_all{n_dset, n_flip}(:,:,n_ctx) = mean(temp_resp,3);
            end
        end
    end
end

resp_all_pool = cell(1,num_flip);
for n_flip = 1:num_flip
    resp_all_pool{1,n_flip} = cat(1,resp_all{:,n_flip});
end

if add_combined
    resp_all_pool{1,num_flip+1} = cat(1,resp_all_pool{:});
end

num_flip = size(resp_all_pool,2);

num_cells = zeros(1, num_flip);
resp_mean = cell(1, num_flip);
resp_sem = cell(1, num_flip);
y_lim_max = 0;
y_lim_min = 0;
for n_flip = 1:num_flip
    num_cells(n_flip) = size(resp_all_pool{n_flip},1);
    resp_mean{n_flip} = squeeze(mean(resp_all_pool{n_flip},1));
    resp_sem{n_flip} = squeeze(std(resp_all_pool{n_flip},[],1)/sqrt(max(num_cells(n_flip)-1,1)));
    
    max_vals = resp_mean{n_flip} + resp_sem{n_flip};
    min_vals = resp_mean{n_flip} - resp_sem{n_flip};
    y_lim_max = max([y_lim_max max(max_vals(:))]);
    y_lim_min = min([y_lim_min min(min_vals(:))]);
end

figure;
for n_flip = 1:num_flip
    subplot(num_flip,1,n_flip); hold on; axis tight; ylim([y_lim_min, y_lim_max]);
    if num_cells(n_flip)
        for n_ctx = 1:size(ctx_plot_list,1)
            color2 = app.ops.context_types_all_colors2{ctx_plot_list(n_ctx,1)};
            shadedErrorBar_YS(plot_t, resp_mean{n_flip}(:,n_ctx),resp_sem{n_flip}(:,n_ctx), color2);
        end

    end

%     if rem(n_ctx,n) ~= 1
%         set(gca,'ytick',[])
%     end
%     if rem(n_ctx,n) == 1
%         if n_ctx > n
%             if app.ConverttoZCheckBox.Value
%                 ylabel(sprintf('%s; %s',app.ops.context_types_labels{mmn_freq(1)}, 'Z scores'));
%             else
%                 ylabel(sprintf('%s; %s',app.ops.context_types_labels{mmn_freq(1)}, 'norm resp'));
%             end
%         else
%             if app.ConverttoZCheckBox.Value
%                 ylabel(sprintf('%s; %s',app.ops.context_types_labels{mmn_freq(2)}, 'Z scores'));
%             else
%                 ylabel(sprintf('%s; %s',app.ops.context_types_labels{mmn_freq(2)}, 'norm resp'));
%             end
%         end
%     end
    %title(sprintf('%s', title2))
end
sgtitle(title_tag, 'Interpreter', 'none')


end