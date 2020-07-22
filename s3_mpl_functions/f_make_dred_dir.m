function dr_params = f_make_dred_dir(dr_params, ops)

%% dim red data saving

% 
% if ops.dred_params.dred_mmn == 1
%     trial_type_tag = 'mmn1_C_fR_D';       
% elseif ops.dred_params.dred_mmn == 2
%     trial_type_tag = 'mmn2_fC_R_fD';   
% elseif ops.dred_params.dred_mmn == 3
%     trial_type_tag = 'mmn12_C_fR_D_fC_R_fD';
% else
%     tt_to_dred = ops.dred_params.trial_types_to_dred;
%     trial_type_tag = num2str(tt_to_dred(1));
%     for n_tr = 2:(numel(tt_to_dred)-1)
%         if and((tt_to_dred(n_tr+1) - tt_to_dred(n_tr)) > 1,(tt_to_dred(n_tr) - tt_to_dred(n_tr-1)) == 1)
%             trial_type_tag = [trial_type_tag  'to' num2str(tt_to_dred(n_tr))];
%         end
%         if (tt_to_dred(n_tr) - tt_to_dred(n_tr-1)) > 1
%             trial_type_tag = [trial_type_tag '_' num2str(tt_to_dred(n_tr))];
%         end
%     end
%     trial_type_tag = ['trials' trial_type_tag 'to' num2str(tt_to_dred(end))];
% end


if ops.dred_params.use_responsive_cells
    resp_tag = 'resp';
else
    resp_tag = '';
end

dr_params.resp_tag = resp_tag;

% for n_cond = 1:numel(ops.regions_to_analyze)
%     cond_name = ops.regions_to_analyze{n_cond};
%     for n_dset = 1:numel(ops.file_names.(cond_name))
save_dred_dir = sprintf('\\%s\\run%.3d_%s_%s\\%s\\%s%s', ops.dred_params.saved_data_dir,ops.dred_params.run_idx, dr_params.trial_type_tag,resp_tag, dr_params.cond_name,ops.paradigm_type,ops.file_names.(dr_params.cond_name){dr_params.n_dset});
if ~isfolder([ops.file_dir save_dred_dir])           % 
    mkdir([ops.file_dir save_dred_dir])
end
dr_params.save_path = save_dred_dir;
%     end
% end

end