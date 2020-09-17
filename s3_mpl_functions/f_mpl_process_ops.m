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
ops.context_types_all= [1:10, 200+(1:7),150, 260, 170, 100+(1:7),250, 160, 270]';
ops.context_types_labels = cell(30,1);
for ii = 1:10
    ops.context_types_labels{ii} = sprintf('Freq %d', ii);
end
for ii = 1:7
    ops.context_types_labels{ii+10} = sprintf('Redf %d', ii);
end
ops.context_types_labels{18} = 'Cont';
ops.context_types_labels{19} = 'RedF pool';
ops.context_types_labels{20} = 'Dev';
for ii = 1:7
    ops.context_types_labels{ii+20} = sprintf('red %d', ii);
end
ops.context_types_labels{28} = 'ContF';
ops.context_types_labels{29} = 'Red pool';
ops.context_types_labels{30} = 'DevF';

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
%ops.context_colors = {'k', 'b', 'r'};
%ops.context_colors = {'k', 'c', 'm'};
%ops.context_colors = ['k', 'b', 'r'];
ops.cond_colors = {[1 0 1], [1 .6, .2], [.2 .8 .2], [0, .6, 1]};
%ops.context_colors = {[0 0 0], [0 1 0], [1 0 1]};
ops.context_colors = {[0 0 0], [0 0 1], [1 0 0]};

ops.context_types_all_colors = zeros(30,1,3);
ops.context_types_all_colors(1:10,:,:) = jet(10);
%ops.context_types_all_colors(1:10,:,:) = reshape(parula(10),10,1,3);
ops.context_types_all_colors(11:17,:,:) = ops.context_colors{2}.*(linspace(0.4,0.8,7))'*0.8;
%ops.context_types_all_colors(11:18,:,:) = parula(8);
ops.context_types_all_colors(18,:,:) = ops.context_colors{1}*0.8;
ops.context_types_all_colors(19,:,:) = ops.context_colors{2}*0.8;
ops.context_types_all_colors(20,:,:) = ops.context_colors{3}*0.8;
ops.context_types_all_colors(21:27,:,:) = ops.context_colors{2}.*(linspace(0.4,0.8,7))';
ops.context_types_all_colors(28,:,:) = ops.context_colors{1};
ops.context_types_all_colors(29,:,:) = ops.context_colors{2};
ops.context_types_all_colors(30,:,:) = ops.context_colors{3};

ops.context_types_all_colors2 = cell(30,1);
for n_cl = 1:10
    ops.context_types_all_colors2{n_cl} = squeeze(ops.context_types_all_colors(n_cl,:,:));
end
for n_cl = 11:17
    ops.context_types_all_colors2{n_cl} = ops.context_colors{2};
end
ops.context_types_all_colors2{18} = ops.context_colors{1};
ops.context_types_all_colors2{19} = ops.context_colors{2};
ops.context_types_all_colors2{20} = ops.context_colors{3};
for n_cl = 21:27
    ops.context_types_all_colors2{n_cl} = ops.context_colors{2};
end
ops.context_types_all_colors2{28} = ops.context_colors{1};
ops.context_types_all_colors2{29} = ops.context_colors{2};  
ops.context_types_all_colors2{30} = ops.context_colors{3};
%figure; imagesc(ops.context_types_all_colors);

ops.context_name = {'Control', 'Redundant', 'Deviant'};
ops.context_name_full = {'Cont', 'RedF', 'Dev', 'ContF', 'Red', 'DevF'};

ops.colors_list2 = {[0, .447, .741], [.85, .325, .098], [.929, .694, .125]...
                   [.4940, .184, .556], [.466, .674, .188]...
                   [.3010, .745, .933], [.635, .078, .184]};
ops.colors_list2 = repmat(ops.colors_list2,1,10);

ops.colors_list = {[1 0 0], [0 1 0], [0 0 1],...
                   [0 1 1], [1 0 1], [1 1 0],...
                   [0 0 0]};
ops.colors_list = repmat(ops.colors_list,1,10);

%ops.fig_title_exp = {'A1 Tones', 'A1 Frequency Grating', 'Combined'};
% MMN order is control, deviant, redundant
% flipMMN order is control, redundant, deviant 
ops.fig_title_run = {'MMN', 'flipMMN', 'Combined'};
ops.errors = {};

ops.plot_params.color_on = [.32 .69 .91];
ops.plot_params.color_off = [.32 0 .91];
ops.plot_params.color_intersect = [0 0 .91];


ops.plot_params.c_map_cont = jet(10);
ops.plot_params.c_map_mmn = [.5 .5 .5; .5 .5 1; 1 .5 .5; .7 .7 .7; .7 .7 1; 1 .7 .7];

end