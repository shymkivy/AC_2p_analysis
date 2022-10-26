function f_dv_plot_mmn(app)

add_combined = 1;

n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = numel(data.experiment);

trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

num_t = sum(trial_frames);

ctx_plot_list = [18, 19, 20; ...
                     28, 29, 30]';
                 
num_flip = size(ctx_plot_list,2);
                 
params = f_dv_gather_params(app);

resp_all = cell(num_dsets, num_flip);
cell_counts = zeros(num_dsets, num_flip);

if strcmpi(app.regiontoplotDropDown.Value, 'all')
    region_num = 0;
elseif strcmpi(app.regiontoplotDropDown.Value, 'A1')
    region_num = 1;
elseif strcmpi(app.regiontoplotDropDown.Value, 'A2')
    region_num = 2;
elseif strcmpi(app.regiontoplotDropDown.Value, 'AAF')
    region_num = 3;
elseif strcmpi(app.regiontoplotDropDown.Value, 'UF')
    region_num = 4;
end


for n_flip = 1:num_flip
    tn_all = ctx_plot_list(:,n_flip);
    ctx1 = app.ops.context_types_all(tn_all)';
    
    for n_dset = 1:num_dsets
        fprintf('dset %d\n', n_dset);
        data1 =  data(n_dset,:);
        stats1 = data1.stats{n_pl};
        params.n_dset = find(data1.idx == app.data.idx);

        cdata = f_dv_compute_cdata(data1, params);

        firing_rate = cdata.S_sm;
        trial_types = data1.trial_types{1};
        stim_times = data1.stim_frame_index{n_pl};
        mmn_freq = data1.MMN_freq{1};
        
        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn_all);
        %cell_is_resp = stats1.peak_resp_cells(:,tn_all);
        
        if ~region_num
            reg_cell_idx = ones(cdata.num_cells,1);
        else
            if and(app.UseregdatalabelsCheckBox.Value, ~isempty(data1.registered_data{1}))
                reg_cell_idx = data1.registered_data{1}.reg_labels==region_num;
            else
                reg_cell_idx = zeros(cdata.num_cells,1);
            end
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
        
        % get resp cells
        resp_cell_idx = logical(sum(resp_cells,2).*reg_cell_idx);
        cell_counts(n_dset, n_flip) = sum(resp_cell_idx);
        
        resp_all{n_dset, n_flip} = zeros(cell_counts(n_dset, n_flip), num_t, size(ctx1,2));
        for n_ctx = 1:size(ctx1,2)
            ctx2 = ctx1(:,n_ctx)';
            temp_resp = trial_data_sort_wctx(resp_cell_idx,:,logical(sum(trial_types_wctx == ctx2,2)));
            resp_all{n_dset, n_flip}(:,:,n_ctx) = mean(temp_resp,3);
            if sum(sum(isnan(resp_all{n_dset, n_flip}(:,:,n_ctx))))
                1;
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
    resp_mean{n_flip} = squeeze(nanmean(resp_all_pool{n_flip},1));
    resp_sem{n_flip} = squeeze(nanstd(resp_all_pool{n_flip},[],1)/sqrt(max(num_cells(n_flip)-1,1)));
    
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
    if app.ConverttoZCheckBox.Value
        ylabel('Z-score');
    else
        ylabel('response mag');
    end
    xlabel('Time (sec)');
    title(sprintf('%s; %d cells', app.regiontoplotDropDown.Value, num_cells(n_flip)));
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
sgtitle([title_tag '; region ' app.regiontoplotDropDown.Value], 'Interpreter', 'none');


end