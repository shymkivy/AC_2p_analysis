function ops = f_mpl_process_ops(ops)
ops.conditions = fieldnames(ops.file_names);
% check if conditions have any datasets

%% remove empty or anavailable datasets
remove_cond = false(numel(ops.conditions),1);
for n_cond = 1:numel(ops.regions_to_analyze)  
    if ~sum(strcmp(ops.regions_to_analyze{n_cond}, ops.conditions))
        warning(['Region ' ops.regions_to_analyze{n_cond} ' is not found in file_names'])
        remove_cond(n_cond) = 1;
    elseif isempty(ops.file_names.(ops.regions_to_analyze{n_cond}))
        disp(['Region ' ops.regions_to_analyze{n_cond} ' has no data and will be removed from analysis'])
        remove_cond(n_cond) = 1;
    end
end
ops.regions_to_analyze(remove_cond) = [];  

%% compute frequencies used
ops.control_carrier_freq = zeros(1, ops.stim.num_freqs);
ops.control_carrier_freq(1) = ops.stim.start_freq;
for ii = 2:ops.stim.num_freqs
    ops.control_carrier_freq(ii) = ops.control_carrier_freq(ii-1) * ops.stim.increase_factor;
end

%%
% create freq legend for plots
ops.context_type_legend = cell(10,3);
for ii = 1:numel(ops.control_carrier_freq)
    ops.context_type_legend{ii,1} = sprintf('%.1fkHz',ops.control_carrier_freq(ii)/1000);
end
ops.context_type_legend{1,2} = 'Cont';
ops.context_type_legend{1,3} = 'ContFlip';
ops.context_type_legend{10,2} = 'Dev';
ops.context_type_legend{10,3} = 'DevFlip';
for ii = 1:8
    ops.context_type_legend{ii+1,2} = sprintf('Red%d', ii);
    ops.context_type_legend{ii+1,3} = sprintf('Red%dFlip', ii);
end

%%
% added 260 and 160 for red sum
ops.context_types_all= [1:10, 200+(1:8), 260, 170, 100+(1:8), 160, 270]';
ops.context_types_labels = cell(28,1);
for ii = 1:10
    ops.context_types_labels{ii} = sprintf('freq %d', ii);
end
for ii = 1:8
    ops.context_types_labels{ii+10} = sprintf('redf %d', ii);
end
ops.context_types_labels{19} = 'redf pool';
ops.context_types_labels{20} = 'dev';
for ii = 1:8
    ops.context_types_labels{ii+20} = sprintf('red %d', ii);
end
ops.context_types_labels{29} = 'red pool';
ops.context_types_labels{30} = 'devf';

%% dim red data saving

if ops.dred_params.dred_mmn == 1
    trial_type_tag = 'mmn1_C_fR_D';       
elseif ops.dred_params.dred_mmn == 2
    trial_type_tag = 'mmn2_fC_R_fD';   
elseif ops.dred_params.dred_mmn == 3
    trial_type_tag = 'mmn12_C_fR_D_fC_R_fD';
else
    tt_to_dred = ops.dred_params.trial_types_to_dred;
    trial_type_tag = num2str(tt_to_dred(1));
    for n_tr = 2:(numel(tt_to_dred)-1)
        if and((tt_to_dred(n_tr+1) - tt_to_dred(n_tr)) > 1,(tt_to_dred(n_tr) - tt_to_dred(n_tr-1)) == 1)
            trial_type_tag = [trial_type_tag  'to' num2str(tt_to_dred(n_tr))];
        end
        if (tt_to_dred(n_tr) - tt_to_dred(n_tr-1)) > 1
            trial_type_tag = [trial_type_tag '_' num2str(tt_to_dred(n_tr))];
        end
    end
    trial_type_tag = ['trials' trial_type_tag 'to' num2str(tt_to_dred(end))];
end
ops.dred_params.trial_type_tag = trial_type_tag;

if ops.dred_params.use_responsive_cells
    resp_tag = 'resp';
else
    resp_tag = '';
end

ops.dred_params.resp_tag = resp_tag;

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    for n_dset = 1:numel(ops.file_names.(cond_name))
        save_dred_dir = sprintf('\\%s\\run%.3d_%s_%s\\%s\\%s%s', ops.dred_params.saved_data_dir,ops.dred_params.run_idx, trial_type_tag,resp_tag, cond_name,ops.paradigm_type,ops.file_names.(cond_name){n_dset});
        if ~isfolder([ops.file_dir save_dred_dir])           % 
            mkdir([ops.file_dir save_dred_dir])
        end
        ops.dred_params.save_path.(cond_name){n_dset} = save_dred_dir;
    end
end

%% subplots dimensions
if numel(ops.regions_to_analyze) == 1
    sm = 1;
    sn = 1;
elseif numel(ops.regions_to_analyze) == 2
    sm = 1;
    sn = 2;
elseif numel(ops.regions_to_analyze) == 3
    sm = 1;
    sn = 3;
elseif numel(ops.regions_to_analyze) == 4
    sm = 2;
    sn = 2;
end
ops.plot_params.reg_sm = sm;
ops.plot_params.reg_sn = sn;

%% other
% ----------------------------figure legends----------------------------
ops.context_colors = {'k', 'b', 'r'};
%ops.context_colors = {'k', 'c', 'm'};
%ops.context_colors = ['k', 'b', 'r'];
ops.context_colors2 = {[0 0 0], [0 1 0], [1 0 1]};
ops.context_name = {'Control', 'Redundant', 'Deviant'};
ops.context_name_full = {'Cont', 'RedF', 'Dev', 'Cont2', 'Red', 'DevF'};

%ops.fig_title_exp = {'A1 Tones', 'A1 Frequency Grating', 'Combined'};
% MMN order is control, deviant, redundant
% flipMMN order is control, redundant, deviant 
ops.fig_title_run = {'MMN', 'flipMMN', 'Combined'};
ops.errors = {};

ops.plot_params.color_on = [83, 177, 232]/256;
ops.plot_params.color_off = [83, 0, 232]/256;
ops.plot_params.color_intersect = [0, 0, 232]/256;


ops.plot_params.c_map_cont = jet(10);
ops.plot_params.c_map_mmn = [.5 .5 .5; .5 .5 1; 1 .5 .5; .7 .7 .7; .7 .7 1; 1 .7 .7];

end