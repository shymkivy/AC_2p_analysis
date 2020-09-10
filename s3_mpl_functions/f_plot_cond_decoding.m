function f_plot_cond_decoding(dec_data_out, plot_var, params, ops)

marginalize_params = 0; % 0= plot for the middle param

%% select params that vary
params_fields = fields(params);
varying_param_bool = false(numel(params_fields),1);
for n_pr = 1:numel(params_fields)
    if ~(strcmpi(params_fields{n_pr}, 'n_rep') || strcmpi(params_fields{n_pr}, plot_var))
        if isnumeric(params.(params_fields{n_pr})) || iscell(params.(params_fields{n_pr}))
            if numel(params.(params_fields{n_pr})) > 1
                varying_param_bool(n_pr) = 1;
            end
        end
    end
end

%% select value to plot
if ~marginalize_params
    varying_param_val = zeros(numel(params_fields),1);
    for n_pr = 1:numel(params_fields)
        if varying_param_bool(n_pr)
            [~, par_ind] = min(abs(params.(params_fields{n_pr}) - mean(params.(params_fields{n_pr}))));
            varying_param_val(n_pr) = params.(params_fields{n_pr})(par_ind);
        end
    end
end
%% plot
if numel(params.(plot_var))>1
    figure; 
    s1= subplot(2,1,1); hold on; axis tight;
    pl = cell(numel(ops.regions_to_analyze),1);
end
dec_data_means = cell(numel(ops.regions_to_analyze),1);

for n_cond = 1:numel(ops.regions_to_analyze)
    dec_data_means{n_cond} = nan(numel(dec_data_out{n_cond}),numel(params.(plot_var)));
    for n_dset = 1:numel(dec_data_out{n_cond})
        for n_param = 1:numel(params.(plot_var))
            temp_st2 = dec_data_out{n_cond}{n_dset}([dec_data_out{n_cond}{n_dset}.(plot_var)] == params.(plot_var)(n_param));
            if ~marginalize_params
                for n_pr = 1:numel(params_fields)
                    if varying_param_bool(n_pr)
                        temp_st2 = temp_st2([temp_st2.(params_fields{n_pr})] == varying_param_val(n_pr));
                    end
                end
            end
            dec_data_means{n_cond}(n_dset,n_param) = mean([temp_st2.accuracy]);
        end
        if numel(params.(plot_var))>1
            p1 = plot(params.(plot_var), dec_data_means{n_cond}(n_dset,:), 'color', ops.cond_colors{n_cond});
            p1.Color(4) = 0.5;
        end
    end
    if numel(params.(plot_var))>1
        pl{n_cond} = plot(params.(plot_var), nanmean(dec_data_means{n_cond}), 'color', ops.cond_colors{n_cond}, 'LineWidth', 2);
    end
end


%figure; hold on;

if numel(params.(plot_var))>1
    subplot(2,1,2); hold on; axis tight;
    for n_cond = 1:numel(ops.regions_to_analyze)
        means1 = nanmean(dec_data_means{n_cond});
        sem1 = nanstd(dec_data_means{n_cond})/sqrt(size(dec_data_means{n_cond},1)-1);
        shadedErrorBar_YS(params.(plot_var), means1, sem1, ops.cond_colors{n_cond});
    end
    legend([pl{:}], ops.regions_to_analyze);
    subplot(s1)
else
    f_plot_dset_deets(dec_data_means, ops);
end

end