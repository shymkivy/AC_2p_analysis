function f_dv_decoder_tuning(app)
% best strategy is to use SVD to extract cell signals, and then do PCA
% before decoder, with mean subtraction step. can get up to 90% decodign

% cosine SVM with no dred also works well it seems?


disp('Decoder...')

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = numel(data.experiment);

strat = 'svd';    % full, onset, offset, mean, pca, svd
num_comp_keep = 5;

use_dim_red = 1;
pca_var_thresh = 70;

if sum(strcmpi(strat, {'svd', 'pca'}))
    title_tag = sprintf('cell %s %d comp', strat, num_comp_keep);
else
    title_tag = sprintf('cell %s', strat);
end
if use_dim_red
    title_tag = sprintf('%s; dred dec %d%% var', title_tag, pca_var_thresh);
end


trial_window = f_str_to_array(app.analysis_BaserespwinEditField.Value);
onset_window = f_str_to_array(app.stats_OnsetRespwinEditField.Value);
offset_window = f_str_to_array(app.stats_OffsetRespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);
num_t = sum(trial_frames);
onset_frames = and(plot_t >= onset_window(1), plot_t <= onset_window(2));
offset_frames = and(plot_t >= offset_window(1), plot_t <= offset_window(2));

[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(app);
num_reg = size(region_num,1);

dec_data_out = cell(num_reg,1);

% tn = [3 4;...
%       4 5;...
%       5 6;...
%       6 7];
tn = [3 4 5 6 7];
%tn = [18 20; 28 30];

[num_gr, num_tn] = size(tn);

mouse_id = cell(num_gr, num_dsets);
dset_id = zeros(num_gr, num_dsets);

params = f_dv_gather_params(app);


dec_params.decoder_type = {'bayes_cosine'}; % 'svm_gaussian' 'svm_cosine' 'svm_linear' 'bayes_cosine'
dec_params.n_rep = 1:10;
dec_params.num_cells = 5:10:80;
dec_params.kFold = 5;

dec_params_list = f_build_param_list(dec_params, {'num_cells', 'n_rep', 'decoder_type'});
num_param_el = numel(dec_params_list);

resp_all = cell(num_dsets, num_gr, num_reg, num_tn);
cell_counts = zeros(num_dsets, num_gr, num_reg);
trial_counts = zeros(num_dsets, num_gr, num_tn);
fprintf('loading data dset #/%d: ', num_dsets);
for n_dset = 1:num_dsets
    fprintf('..%d', n_dset);
    data1 =  data(n_dset,:);
    stats1 = data1.stats{n_pl};
    params.n_dset = find(data1.idx == app.data.idx);

    cdata = f_dv_compute_cdata(data1, params);

    firing_rate = cdata.S_sm;
    mmn_freq = data1.MMN_freq{1};
    trial_types = data1.trial_types{1};
    trial_types_ctx = f_dv_mark_tt_ctx(trial_types, mmn_freq, app.ops);

    stim_times = data1.stim_frame_index{n_pl};
    
    trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);

    reg_cell_labels = f_dv_get_area_label(app, data1);

    for n_gr = 1:num_gr
        mouse_id{n_gr, n_dset} = data1.mouse_id{1};
        dset_id(n_gr, n_dset) = n_dset;
        
        tn1 = tn(n_gr,:);
        
        resp_cells = f_dv_get_resp_vals_cells(app, stats1, tn1);
        resp_cells2 = sum(resp_cells,2);
        %cell_is_resp = stats1.peak_resp_cells(:,tn_all);

        for n_reg = 1:num_reg
            n_reg2 = region_num(n_reg,:);
            
            % get resp cells
            reg_cell_idx = logical(sum(reg_cell_labels == n_reg2,2));
            resp_cell_idx = logical(resp_cells2.*reg_cell_idx);
            
            num_cells = sum(resp_cell_idx);
            
            cell_counts(n_dset, n_gr, n_reg) = num_cells;
            for n_ctx = 1:num_tn
                ctx1 = app.ops.context_types_all(tn1(n_ctx));
                %idx1 = trial_types_wctx == ctx2;
                idx1 = logical(sum([trial_types, trial_types_ctx] == ctx1,2));
                num_trials = sum(idx1);
                if num_cells
                    if strcmpi(strat, 'onset')
                        resp_all{n_dset, n_gr, n_reg, n_ctx} = reshape(mean(trial_data_sort(resp_cell_idx,onset_frames,idx1),2), num_cells, num_trials);
                    elseif strcmpi(strat, 'offset')
                        resp_all{n_dset, n_gr, n_reg, n_ctx} = reshape(mean(trial_data_sort(resp_cell_idx,offset_frames,idx1),2), num_cells, num_trials);
                    elseif strcmpi(strat, 'mean')
                        resp_all{n_dset, n_gr, n_reg, n_ctx} = reshape(mean(trial_data_sort(resp_cell_idx,:,idx1),2), num_cells, num_trials);
                    elseif strcmpi(strat, 'pca')
                        resp_all{n_dset, n_gr, n_reg, n_ctx} = zeros(num_cells, num_comp_keep, num_trials);
                        temp_data = trial_data_sort(resp_cell_idx,:,idx1);
                        for n_cell = 1:num_cells
                            temp_data2 = squeeze(temp_data(n_cell,:,:))';
                            [~,score,~,~,explained,mu] = pca(temp_data2);
                            %disp(explained(1))
                            resp_all{n_dset, n_gr, n_reg, n_ctx}(n_cell,:,:) = score(:,1:num_comp_keep)';
                        end
                    elseif strcmpi(strat, 'svd')
                        resp_all{n_dset, n_gr, n_reg, n_ctx} = zeros(num_cells, num_comp_keep, num_trials);
                        temp_data = trial_data_sort(resp_cell_idx,:,idx1);
                        for n_cell = 1:num_cells
                            temp_data2 = squeeze(temp_data(n_cell,:,:))';
                            %temp_data2 = temp_data2 - mean(temp_data2);
                            [U, S, V] = svd(temp_data2,"econ");
                            explained = diag(S).^2 / sum(diag(S).^2)*100;
                            %disp(explained(1))
                            score = U*S;
                            resp_all{n_dset, n_gr, n_reg, n_ctx}(n_cell,:,:) = score(:,1:num_comp_keep)';
                        end
                    elseif strcmpi(strat, 'full')
                        resp_all{n_dset, n_gr, n_reg, n_ctx} = trial_data_sort(resp_cell_idx,:,idx1);
                    end
                end
                trial_counts(n_dset, n_gr, n_ctx) = num_trials;
            end
        end
    end
end
fprintf('\n')

num_dims = ndims(resp_all{1, 1, 1,1});

fprintf('decoder dset #/%d: ', num_dsets);
dec_data1 = cell(num_dsets,num_gr,num_reg);
for n_dset = 1:num_dsets
    fprintf('..%d', n_dset);
    for n_reg = 1:num_reg
    reg_name = app.ops.regions_to_analyze{n_reg};
        for n_tt_gr = 1:num_gr
            trial_counts2 = squeeze(trial_counts(n_dset, n_tt_gr,:));
            trials_use = min(trial_counts2);

            tn1 = tn(n_tt_gr,:);
            tt1 = app.ops.context_types_all(tn1)';
            
            trial_types1 = reshape(ones(trials_use, num_tn) .* tt1 , trials_use*num_tn,1);

            temp_params = dec_params_list;
            for n_param_el = 1:num_param_el
                temp_params(n_param_el).cond_name = reg_name;
                temp_params(n_param_el).n_dset = n_dset;
                temp_params(n_param_el).tt = tt1;

                traces0 = cell(num_tn, 1);
                for n_tn = 1:num_tn
                    samp_tr = randsample(trial_counts2(n_tn), trials_use, false);
                    traces00 = resp_all{n_dset, n_tt_gr, n_reg, n_tn};
                    if ~isempty(traces00)
                        if num_dims == 3
                            traces0{n_tn} = traces00(:,:,samp_tr);
                        elseif num_dims == 2
                            traces0{n_tn} = traces00(:,samp_tr);
                        end
                    end
                end
                
                if num_dims == 3
                    traces1 = cat(3,traces0{:});
                    [num_cells, num_t, num_trials] = size(traces1);
                elseif num_dims == 2
                    traces1 = cat(2,traces0{:});
                    [num_cells, num_trials] = size(traces1);
                end
                
                if temp_params(n_param_el).num_cells < num_cells
    
                    % randomize trial order
                    rand_tr_order = randperm(num_trials);
                    rand_tr_order_shuff = randperm(num_trials);
                    trial_types2 = trial_types1(rand_tr_order);
                    trial_types_shuff = trial_types2(rand_tr_order_shuff);
        
                    if num_dims == 3
                        traces2 = traces1(:,:,rand_tr_order);
                    elseif num_dims == 2
                        traces2 = traces1(:,rand_tr_order);
                    end
        
                    %
                    response = trial_types2;
                    response_shuff = trial_types_shuff;

                    cells_pred = randsample(num_cells, temp_params(n_param_el).num_cells);

                    if num_dims == 3
                        traces3 = reshape(traces2(cells_pred,:,:), temp_params(n_param_el).num_cells*num_t, []);
                    elseif num_dims == 2
                        traces3 = traces2(cells_pred,:);
                    end

                    if use_dim_red
                    % reduce dim with PCA
                        % seems like removing mean with PCA improves performance
                        %[U, S, V] = svd(traces3',"econ");
                        %score = U*S;
                        %explained = diag(S).^2 / sum(diag(S).^2)*100;
                        [~,score,~,~,explained,mu] = pca(traces3');
                        num_comps = max(sum(cumsum(explained)<pca_var_thresh),1);
                        predictors = score(:,1:num_comps);
                    else
                        predictors = traces3';
                    end
                    %%
                    if strcmpi(temp_params(n_param_el).decoder_type, 'svm_gaussian')
                        temp_params(n_param_el).KernelFunction = 'gaussian';
                        temp_params(n_param_el).KernelScale = 5.5;
                        temp_params(n_param_el).accuracy = f_svm_decoder(predictors, response, tt1, temp_params(n_param_el));
                        temp_params(n_param_el).accuracy_shuff = f_svm_decoder(predictors, response_shuff, tt1, temp_params(n_param_el));
                    elseif strcmpi(temp_params(n_param_el).decoder_type, 'svm_cosine')
                        temp_params(n_param_el).KernelFunction = 'cosineKernel';
                        temp_params(n_param_el).accuracy = f_svm_decoder(predictors, response, tt1, temp_params(n_param_el));
                        temp_params(n_param_el).accuracy_shuff = f_svm_decoder(predictors, response_shuff, tt1, temp_params(n_param_el));
                    elseif strcmpi(temp_params(n_param_el).decoder_type, 'svm_linear')
                        temp_params(n_param_el).KernelFunction = 'linear';
                        temp_params(n_param_el).accuracy = f_svm_decoder(predictors, response, tt1, temp_params(n_param_el));
                        temp_params(n_param_el).accuracy_shuff = f_svm_decoder(predictors, response_shuff, tt1, temp_params(n_param_el));
                    elseif strcmpi(temp_params(n_param_el).decoder_type, 'bayes_cosine')
                        temp_params(n_param_el).KernelFunction = 'cosine';
                        temp_params(n_param_el).accuracy = f_bayes_decoder_wrap(predictors, response, tt1, temp_params(n_param_el));
                        temp_params(n_param_el).accuracy_shuff = f_bayes_decoder_wrap(predictors, response_shuff, tt1, temp_params(n_param_el));
                    end
                else
                    temp_params(n_param_el).accuracy = NaN;
                    temp_params(n_param_el).accuracy_shuff = NaN;
                end
            end
            %dec_data1{n_dset, n_tt_cat} = temp_params;
            dec_data1{n_dset, n_tt_gr, n_reg} = temp_params;
        end
    end
end
fprintf('\n')
for n_reg = 1:num_reg
    dec_data2 = dec_data1(:,:,n_reg);
    dec_data_out{n_reg} = dec_data2(:);
end
%
%% plot crap

if numel(dec_params.num_cells) > 1
    f_plot_cond_decoding(dec_data_out, 'num_cells', dec_params, app.ops, title_tag);
    if numel(dec_params.decoder_type) > 1
        f_plot_decoret_types(dec_data_out, 'num_cells', dec_params, app.ops, title_tag);
    end
end

if numel(dec_params.kFold) > 1
    f_plot_cond_decoding(dec_data_out, 'kFold', dec_params, app.ops);
end

disp('Done')


end