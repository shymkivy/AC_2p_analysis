function dec_data = f_decoder_binwise_onevall(data_all, tt_all, params)

num_dsets = numel(data_all);
[~, num_frames, ~] = size(data_all{1});
num_trials = numel(unique(tt_all{1}));

fprintf('dset n/%d: ', num_dsets);

dec_acc_frames = zeros(num_dsets, num_frames);
dec_acc_frames_shuff = zeros(num_dsets, num_frames);
dec_acc_frames_bycl = zeros(num_dsets, num_frames, num_trials);
dec_acc_frames_bycl_shuff = zeros(num_dsets, num_frames, num_trials);

for n_dset = 1:num_dsets
    
    fprintf('%d..', n_dset);
    
    trial_data_sort3 = data_all{n_dset};
    trial_types3 = tt_all{n_dset};
                
    sh_idx1 = randperm(numel(trial_types3));
    
    trial_types4 = trial_types3(sh_idx1);
    trial_data_sort4 = trial_data_sort3(:,:,sh_idx1);

    trial_types4_shuff = trial_types4(randperm(numel(trial_types4)));
    
    for n_fr = 1:sum(num_frames)
        
        trial_data_sort4n = squeeze(trial_data_sort4(:,n_fr,:));
        trial_data_sort4n = trial_data_sort4n - mean(trial_data_sort4n(:));
        trial_data_sort4n = trial_data_sort4n/std(trial_data_sort4n(:));
    
        if strcmpi(params.decoder_type, 'tree')
            dec_out = decoder_tree(trial_data_sort4n', trial_types4);
            dec_out_shuff = decoder_tree(trial_data_sort4n', trial_types4_shuff);
        elseif strcmpi(params.decoder_type, 'bayes')
            dec_out = decoder_naivebayes(trial_data_sort4n', trial_types4);
            dec_out_shuff = decoder_naivebayes(trial_data_sort4n', trial_types4_shuff);
        elseif strcmpi(params.decoder_type, 'svm')
            dec_out = decoder_svm(trial_data_sort4n', trial_types4, 1);
            dec_out_shuff = decoder_svm(trial_data_sort4n', trial_types4_shuff, 1);
        end

        dec_acc_frames(n_dset,n_fr) = dec_out.validationAccuracy;
        dec_acc_frames_bycl(n_dset,n_fr,:) = dec_out.acc_by_class;
        dec_acc_frames_shuff(n_dset,n_fr) = dec_out_shuff.validationAccuracy;
        dec_acc_frames_bycl_shuff(n_dset,n_fr,:) = dec_out_shuff.acc_by_class;
    end

end
dec_data = struct();
dec_data.accuracy = dec_acc_frames;
dec_data.accuracy_by_class = dec_acc_frames_bycl;
dec_data.accuracy_shuff = dec_acc_frames_shuff;
dec_data.accuracy_by_class_shuff = dec_acc_frames_bycl_shuff;

end