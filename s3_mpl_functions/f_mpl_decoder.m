function f_mpl_decoder(data, ops)
disp('Decoder...')
f = waitbar(0,'decoding');

dec_data_out = cell(numel(ops.regions_to_analyze),1);

tt = [3 4 5 6 170 270];
dec_cell_start = 5;
dec_cell_int = 3;
dec_cell_max = 20;
num_reps = 10;

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    dec_data1 = cell(cdata.num_dsets,1);
    for n_dset = 1:cdata.num_dsets
        
        trial_types = cdata.trial_types_pr{n_dset}(1:800);
        traces = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag(:,1:800);
        resp_cells = logical(sum(cdata.peak_tuned_trials_full{n_dset}(:,logical(sum(ops.context_types_all == tt,2))),2));
        traces = traces(resp_cells,:);
        
        tr_ind = logical(sum(trial_types == tt,2));
        traces2 = traces(:,tr_ind);
        
        trial_types2 = trial_types(tr_ind);
        num_cells = size(traces,1);
        max_cells1 = min(dec_cell_max, num_cells);

        response = trial_types2;
        % cycle through cell nums

        cell_list = dec_cell_start:dec_cell_int:max_cells1;
        

        accuracy1 = zeros(numel(cell_list),num_reps);

        for n_celln = 1:numel(cell_list)
            for n_rep = 1:num_reps

                cells1 = randsample(num_cells, cell_list(n_celln));


                [~,score,~,~,explained,~] = pca(traces2(cells1,:)');
                num_comps = max(sum(cumsum(explained)<70),1);
                predictors = score(:,1:num_comps);
                %%
        %         SVMModel = fitcsvm(...
        %             predictors, ...
        %             response, ...
        %             'KernelFunction', 'gaussian', ...     % 'KernelFunction', 'polynomial' 'gaussian'
        %             'PolynomialOrder', [], ...               % 2
        %             'KernelScale', 7.1, ...              % auto fine = 1.8 medium= 7.1 coarse = 28
        %             'BoxConstraint', 1, ...
        %             'Standardize', true, ...
        %             'ClassNames', tt');

                template = templateSVM(...
                    'KernelFunction', 'gaussian', ...
                    'PolynomialOrder', [], ...
                    'KernelScale', 10, ...
                    'BoxConstraint', 1, ...
                    'Standardize', true);
                SVMModel = fitcecoc(...
                    predictors, ...
                    response, ...
                    'Learners', template, ...
                    'Coding', 'onevsone', ...
                    'ClassNames', tt');

                %% cross validation 
                model1 = SVMModel;

                partitionedModel = crossval(model1, 'KFold', 5);

                %[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

                % Compute validation accuracy
                accuracy1(n_celln, n_rep) = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

            end
        end
        dec_data1{n_dset} = accuracy1;
        waitbar(n_dset/cdata.num_dsets,f,cond_name);
    end
    dec_data_out{n_cond} = dec_data1;
    
%     figure(f1);
%     shadedErrorBar(cell_list, mean(accuracy1,2), std(accuracy1,[],2)./sqrt(num_reps-1))
end

close(f);
%% plot crap
params_dec.dec_cell_start = dec_cell_start;
params_dec.dec_cell_int = dec_cell_int;
params_dec.dec_cell_max = dec_cell_max;
params_dec.num_reps = num_reps;

f_plot_cont_decoding(dec_data_out, params_dec, ops)


disp('Done')

end