function [tn_to_dred, trial_type_tag] = f_select_trial_type(trial_types_to_dred, cdata,n_dset, ops)

if isnumeric(trial_types_to_dred)
    tn_to_dred = trial_types_to_dred;
elseif strcmpi(trial_types_to_dred, 'all')
    tn_to_dred = 1:numel(ops.context_types_all);
elseif strcmpi(trial_types_to_dred, 'mmn1')
    tn_to_dred = cdata.ctx_mmn{n_dset}(1:3);
elseif strcmpi(trial_types_to_dred, 'mmn2')
    tn_to_dred = cdata.ctx_mmn{n_dset}(4:6);
elseif strcmpi(trial_types_to_dred, 'mmn12')
    tn_to_dred = cdata.ctx_mmn{n_dset};
elseif strcmpi(trial_types_to_dred, 'cont_all')
    tn_to_dred = 1:10;
elseif strcmpi(trial_types_to_dred, 'dd1')
    tn_to_dred = cdata.ctx_mmn{n_dset}(3);
elseif strcmpi(trial_types_to_dred, 'dd2')
    tn_to_dred = cdata.ctx_mmn{n_dset}(6);
elseif strcmpi(trial_types_to_dred, 'dd12')
    tn_to_dred = cdata.ctx_mmn{n_dset}([3,6]);
elseif strcmpi(trial_types_to_dred, 'red1')
    tn_to_dred = cdata.ctx_mmn{n_dset}(2);
elseif strcmpi(trial_types_to_dred, 'red2')
    tn_to_dred = cdata.ctx_mmn{n_dset}(5);
elseif strcmpi(trial_types_to_dred, 'red12')
    tn_to_dred = cdata.ctx_mmn{n_dset}([2,5]);
elseif strcmpi(trial_types_to_dred, 'cont1')
    tn_to_dred = cdata.ctx_mmn{n_dset}(1);
elseif strcmpi(trial_types_to_dred, 'cont2')
    tn_to_dred = cdata.ctx_mmn{n_dset}(4);
elseif strcmpi(trial_types_to_dred, 'cont12')
    tn_to_dred = cdata.ctx_mmn{n_dset}([1,4]);
end

if isnumeric(trial_types_to_dred)
    trial_type_tag = ['trials ' num2str(trial_types_to_dred)];
else
    trial_type_tag = trial_types_to_dred;
end


end