function dred_data_list = f_dim_red_cv(trial_data, ops, cv_params)
cond_name = cv_params.cond_name;
n_dset = cv_params.n_dset;
volume_period = cv_params.volume_period;
tt_to_dred = cv_params.tt_to_dred;

[num_cells, ~, num_trials] = size(trial_data);


%% make sure every fold has active cells, need positive definite matrix
active_cells = true(num_cells,1);
num_folds = ops.dred_params.cv_num_folds;
cv_bins = floor(linspace(1,num_trials+1,num_folds+1));
% find active cells 
test_data_ind = false(num_trials,num_folds);
for n_cv = 1:num_folds
    test_data_ind(cv_bins(n_cv):(cv_bins(n_cv+1)-1),n_cv) = 1;
    active_cells = active_cells.*sum(sum(trial_data(:,:,~test_data_ind(:,n_cv)),3),2)>0;
end

trial_data1 = trial_data(active_cells,:,:);

%% create list of data to process
dred_data_list = struct();
dred_idx = 1;
for n_met = 1:numel(ops.dred_params.method_list)
    if strcmpi(ops.dred_params.method_list{n_met}, 'gpfa')
        kern1 = 0;
    else
        kern1 = ops.dred_params.kernSD;
    end
    for n_kern = 1:numel(kern1)
        for n_comp = 1:numel(ops.dred_params.num_comp)
            for n_cv = 1:ops.dred_params.cv_num_folds
                if num_cells>=ops.dred_params.num_comp(n_comp)
                    dred_data_list(dred_idx).method = ops.dred_params.method_list{n_met};
                    dred_data_list(dred_idx).kernSD = kern1(n_kern);
                    dred_data_list(dred_idx).n_comp = ops.dred_params.num_comp(n_comp);
                    dred_data_list(dred_idx).n_cv = n_cv;
                    dred_data_list(dred_idx).cv_num_folds = ops.dred_params.cv_num_folds;
                    dred_data_list(dred_idx).file_name = [cond_name '_' ops.paradigm_type ops.file_names.(cond_name){n_dset}];
                    dred_data_list(dred_idx).vol_period = volume_period;
                    dred_data_list(dred_idx).trial_types = tt_to_dred;
                    dred_data_list(dred_idx).num_cells = num_cells;
                    dred_data_list(dred_idx).cond_name = cond_name;
                    dred_data_list(dred_idx).n_dset = n_dset;
                    dred_idx = dred_idx + 1;
                end
            end
        end
    end
end

%%
save_path = [ops.file_dir cv_params.save_path '\'];
file_list = dir([save_path '*.mat']);
file_list_names = {file_list.name};



for n_dt = 1:numel(dred_data_list)
    dred_save_name = sprintf('%s_s%d_c%d_cv%dof%d', dred_data_list(n_dt).method,...
        dred_data_list(n_dt).kernSD, dred_data_list(n_dt).n_comp,...
        dred_data_list(n_dt).n_cv, dred_data_list(n_dt).cv_num_folds);

    fprintf('%s dim red train file %d of %d %s\n',dred_data_list(n_dt).file_name, n_dt, numel(dred_data_list), dred_save_name);

    yTrain_3d = trial_data1(:,:,~test_data_ind(:,dred_data_list(n_dt).n_cv));
    yTest_3d = trial_data1(:,:,test_data_ind(:,dred_data_list(n_dt).n_cv));
    if ~strcmpi(dred_data_list(n_dt).method,'gpfa')
        yTrain_3ds = f_smooth_gauss(yTrain_3d, dred_data_list(n_dt).kernSD/volume_period);
        yTest_3ds = f_smooth_gauss(yTest_3d, dred_data_list(n_dt).kernSD/volume_period);
    else
        yTrain_3ds = yTrain_3d;
        yTest_3ds = yTest_3d;
    end

    if isempty(file_list_names) || ~sum(strcmpi([dred_save_name '.mat'], file_list_names))
        [dred_factors, ydred_data] = f_dred_train(yTrain_3ds, dred_data_list(n_dt).n_comp, dred_data_list(n_dt).method);
        dred_factors.train_err = norm(yTrain_3d(:) - ydred_data(:))/norm(yTrain_3d(:));
        dred_params = dred_data_list(n_dt);
        save([save_path dred_save_name], 'dred_factors', 'dred_params');
    else
        load([save_path dred_save_name]);
    end

    if ~isfield(dred_factors, 'test_err') || ~isfield(dred_factors, 'test_err_sm')
        Ycs = f_dred_test(yTest_3ds, dred_factors.dred_factors, dred_data_list(n_dt).method);
        dred_factors.test_err = norm(yTest_3d(:) - Ycs(:))/norm(yTest_3d(:));
        dred_factors.test_err_sm = norm(yTest_3ds(:) - Ycs(:))/norm(yTest_3ds(:));
        save([save_path dred_save_name], 'dred_factors', 'dred_params');
    end
    dred_data_list(n_dt).train_err = dred_factors.train_err;
    dred_data_list(n_dt).train_err_sm = dred_factors.train_err_sm;
    dred_data_list(n_dt).test_err = dred_factors.test_err;
    dred_data_list(n_dt).test_err_sm = dred_factors.test_err_sm;
    dred_data_list(n_dt).dred_factors = dred_factors.dred_factors;
end

fig1 = f_plot_dred(dred_data_list);
for n_fig = 1:numel(fig1)
    figure(fig1{n_fig})
    suptitle([cond_name 'dset ' num2str(n_dset) ' n = ' num2str(sum(active_cells))]);
end
end