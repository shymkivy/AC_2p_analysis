function f_mpl_decoder(data, ops)
disp('Decoder...')
f = waitbar(0,'decoding');

dec_data_out = cell(numel(ops.regions_to_analyze),1);

% tn = [3 4;...
%       4 5;...
%       5 6;...
%       6 7];
%tn = [3 4 5 6 7];
tn = [18 20; 28 30];
tt = ops.context_types_all(tn);


random_sample = 0; % 0 = sort and sequentially take
use_dim_red = 1;
sort_mag = 0; % 0 = reliability

dec_params.decoder_type = {'bayes_cosine'}; % 'svm_gaussian' 'svm_cosine' 'svm_linear' 'bayes_cosine'
dec_params.n_rep = 1:10;
dec_params.dec_num_cells = 5:10:80;
dec_params.kFold = 5;


num_tt = numel(tt);

dec_params_list = f_build_param_list(dec_params, {'dec_num_cells', 'n_rep', 'decoder_type'});
num_param_el = numel(dec_params_list);

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data(strcmpi(data.area, cond_name),:);
    
    num_dsets = numel(cdata.area);
    num_tt_cat = size(tt,1);
    dec_data1 = cell(num_dsets,num_tt_cat);
    for n_dset = 1:num_dsets
        for n_tt_cat = 1:num_tt_cat
            trial_types = cdata.trial_types_wctx{n_dset};
            traces = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
            tt2 = tt(n_tt_cat,:);

            if sort_mag
                traces_ave = cdata.trial_ave_z{n_dset};
                resp_mag = squeeze(max(traces_ave, [], 2));
            else
                resp_mag = cdata.peak_tuned_trials_full_reliab{n_dset};
            end

            % select trials
            tr_ind = logical(sum(trial_types == tt2,2));
            traces2 = traces(:,tr_ind);
            trial_types2 = trial_types(tr_ind);

            %
            [num_cells, num_trials] = size(traces2);

            % randomize trial order
            rand_tr_order = randperm(num_trials);
            traces2 = traces2(:,rand_tr_order);
            trial_types2 = trial_types2(rand_tr_order);

            %
            resp_mag2 = resp_mag(:,tn);

            %
            response = trial_types2;
            temp_params = dec_params_list;
            for n_param_el = 1:num_param_el
                temp_params(n_param_el).cond_name = cond_name;
                temp_params(n_param_el).n_dset = n_dset;
                temp_params(n_param_el).tt = tt2;
                if temp_params(n_param_el).dec_num_cells < num_cells
                    if random_sample
                        % choose cells
                        cells_pred = randsample(num_cells, temp_params(n_param_el).dec_num_cells);
                    else
                        list_samp = [(1:num_cells)', resp_mag2];
                        cells_pred = zeros(temp_params(n_param_el).dec_num_cells,1);
                        for n_cell_ind = 1:temp_params(n_param_el).dec_num_cells
                            curr_tt = rem(n_cell_ind-1,num_tt)+1;
                            [~, n_cell_ind2]= max(list_samp(:,curr_tt+1));
                            n_cell = list_samp(n_cell_ind2,1);
                            list_samp(n_cell_ind2,:) = [];
                            cells_pred(n_cell_ind) = n_cell;
                        end
                    end

                    %X1 = [traces2(cells_pred,:)', response]

                    if use_dim_red
                    % reduce dim with PCA
                        [~,score,~,~,explained,~] = pca(traces2(cells_pred,:)');
                        num_comps = max(sum(cumsum(explained)<70),1);
                        predictors = score(:,1:num_comps);
                    else
                        predictors = traces2(cells_pred,:)';
                    end
                    %%
                    if strcmpi(temp_params(n_param_el).decoder_type, 'svm_gaussian')
                        temp_params(n_param_el).KernelFunction = 'gaussian';
                        temp_params(n_param_el).KernelScale = 5.5;
                        temp_params(n_param_el).accuracy = f_svm_decoder(predictors, response, tt2, temp_params(n_param_el));
                    elseif strcmpi(temp_params(n_param_el).decoder_type, 'svm_cosine')
                        temp_params(n_param_el).KernelFunction = 'cosineKernel';
                        temp_params(n_param_el).accuracy = f_svm_decoder(predictors, response, tt2, temp_params(n_param_el));
                    elseif strcmpi(temp_params(n_param_el).decoder_type, 'svm_linear')
                        temp_params(n_param_el).KernelFunction = 'linear';
                        temp_params(n_param_el).accuracy = f_svm_decoder(predictors, response, tt2, temp_params(n_param_el));
                    elseif strcmpi(temp_params(n_param_el).decoder_type, 'bayes_cosine')
                        temp_params(n_param_el).KernelFunction = 'cosine';
                        temp_params(n_param_el).accuracy = f_bayes_decoder_wrap(predictors, response, tt2, temp_params(n_param_el));
                    end
                else
                    temp_params(n_param_el).accuracy = NaN;
                end
            end
            dec_data1{n_dset, n_tt_cat} = temp_params;
        end
        
        waitbar(n_dset/num_dsets,f,['Decoder, ' cond_name]);
    end
    dec_data_out{n_cond} = dec_data1(:);
    
%     figure(f1);
%     shadedErrorBar(cell_list, mean(accuracy1,2), std(accuracy1,[],2)./sqrt(num_reps-1))
end

close(f);
%% plot crap

if numel(dec_params.dec_num_cells) > 1
    f_plot_cond_decoding(dec_data_out, 'dec_num_cells', dec_params, ops);
    if numel(dec_params.decoder_type) > 1
        f_plot_decoret_types(dec_data_out, 'dec_num_cells', dec_params, ops);
    end
end

if numel(dec_params.kFold) > 1
    f_plot_cond_decoding(dec_data_out, 'kFold', dec_params, ops);
end

disp('Done')

end