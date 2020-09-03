function f_mpl_decoder(data, ops)
disp('Decoder...')
f = waitbar(0,'decoding');

dec_data_out = cell(numel(ops.regions_to_analyze),1);

tn = [18 20 28 30];
tt = ops.context_types_all(tn)';
dec_num_cells = 70;

num_reps = 10;

random_sample = 0; % 0 = sort and sequentially take
sort_mag = 1; % 0 = reliability
use_dim_red = 1;

dec_params.KernelFunction = 'gaussian';
dec_params.KernelScale = 7;
dec_params.kFold = 5;

num_tt = numel(tt);

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    dec_data1 = cell(cdata.num_dsets,1);
    for n_dset = 1:cdata.num_dsets
        
        trial_types = cdata.trial_types_wctx{n_dset};
        traces = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
        
        if sort_mag
            traces_ave = cdata.trial_ave_z{n_dset};
            resp_mag = squeeze(max(traces_ave, [], 2));
        else
            resp_mag = cdata.peak_tuned_trials_full_reliab{n_dset};
        end
        
        
        % select trials
        tr_ind = logical(sum(trial_types == tt,2));
        traces2 = traces(:,tr_ind);
        trial_types2 = trial_types(tr_ind);
        resp_mag2 = resp_mag(:,tn);
        
        
%         figure; hold on;
%         for n_c = 1:4
%             ecdf(resp_mag2(:,n_c));
%         end

%         % select cells
%         resp_cells = logical(sum(cdata.peak_tuned_trials_full{n_dset}(:,logical(sum(ops.context_types_all == tt,2))),2));
%         traces3 = traces2(resp_cells,:);
        

        % check if there are enough cells in dset 
        num_cells = size(traces2,1);
        dec_num_cells2 = dec_num_cells(dec_num_cells<num_cells);
        
        response = trial_types2;
        
        accuracy1 = zeros(numel(dec_num_cells2),num_reps);

        for n_celln = 1:numel(dec_num_cells2)
            for n_rep = 1:num_reps
                if random_sample
                    % choose cells
                    cells_pred = randsample(num_cells, dec_num_cells2(n_celln));
                else
                    list_samp = [(1:num_cells)', resp_mag2];
                    cells_pred = zeros(dec_num_cells2(n_celln),1);
                    for n_cell_ind = 1:dec_num_cells2(n_celln)
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
                    num_comps = max(sum(cumsum(explained)<80),1);
                    predictors = score(:,1:num_comps);
                else
                    predictors = traces2(cells_pred,:)';
                end
                %%
                
                accuracy1(n_celln, n_rep) = f_svm_decoder(predictors, response, tt, dec_params);              
            end
        end
        dec_data1{n_dset} = accuracy1;
        waitbar(n_dset/cdata.num_dsets,f,['Decoder, ' cond_name]);
    end
    dec_data_out{n_cond} = dec_data1;
    
%     figure(f1);
%     shadedErrorBar(cell_list, mean(accuracy1,2), std(accuracy1,[],2)./sqrt(num_reps-1))
end

close(f);
%% plot crap
params_dec.dec_num_cells = dec_num_cells;
params_dec.num_reps = num_reps;

f_plot_cond_decoding(dec_data_out, params_dec, ops)
title(['Decoder, conditions ' num2str(tt)])

disp('Done')

end