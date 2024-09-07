function f_plot_cond_decoding(dec_data_out, plot_var, params, ops, title_tag)

marginalize_params = 0; % 0= plot for the middle param

%% select params that vary
params_fields = fields(params);
varying_param_bool = false(numel(params_fields),1);
for n_pr = 1:numel(params_fields)
    if ~(strcmpi(params_fields{n_pr}, 'n_rep') || strcmpi(params_fields{n_pr}, plot_var) || strcmpi(params_fields{n_pr}, 'decoder_type'))
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
for n_dec = 1:numel(params.decoder_type)
    decoder1 = params.decoder_type{n_dec};
    title_tag2 = sprintf('Decoder %s; %s', decoder1, title_tag);
    if numel(params.(plot_var))>1
        figure; 
        s1= subplot(2,1,1); hold on; axis tight;
        pl = cell(numel(ops.regions_to_analyze),1);
    end
    dec_data_means = cell(numel(ops.regions_to_analyze),1);
    dec_data_means_shuff = cell(numel(ops.regions_to_analyze),1);
    for n_reg = 1:numel(ops.regions_to_analyze)
        dec_data_means{n_reg} = nan(numel(dec_data_out{n_reg}),numel(params.(plot_var)));
        dec_data_means_shuff{n_reg} = nan(numel(dec_data_out{n_reg}),numel(params.(plot_var)));
        for n_dset = 1:numel(dec_data_out{n_reg})
            for n_param = 1:numel(params.(plot_var))
                temp_st2 = dec_data_out{n_reg}{n_dset}([dec_data_out{n_reg}{n_dset}.(plot_var)] == params.(plot_var)(n_param));
                temp_st3 = temp_st2(strcmpi([temp_st2.decoder_type],decoder1));
                if ~marginalize_params
                    for n_pr = 1:numel(params_fields)
                        if varying_param_bool(n_pr)
                            temp_st3 = temp_st3([temp_st3.(params_fields{n_pr})] == varying_param_val(n_pr));
                        end
                    end
                end
                dec_data_means{n_reg}(n_dset,n_param) = mean([temp_st3.accuracy]);
                dec_data_means_shuff{n_reg}(n_dset,n_param) = mean([temp_st3.accuracy_shuff]);
            end
            if numel(params.(plot_var))>1
                p2 = plot(params.(plot_var), dec_data_means_shuff{n_reg}(n_dset,:), color='black');
                p2.Color(4) = 0.5;
                p1 = plot(params.(plot_var), dec_data_means{n_reg}(n_dset,:), color=ops.cond_colors{n_reg});
                p1.Color(4) = 0.5;
            end
        end
        if numel(params.(plot_var))>1
            pls = plot(params.(plot_var), mean(dec_data_means_shuff{n_reg},'omitnan'), color='black', LineWidth=2);
            pl{n_reg} = plot(params.(plot_var), mean(dec_data_means{n_reg},'omitnan'), color=ops.cond_colors{n_reg}, LineWidth=2);
        end
    end
    ylabel('performance')
    xlabel(strrep(plot_var,'_',' '))
    title(strrep(title_tag2,'_',' '))

    %figure; hold on;

    if numel(params.(plot_var))>1
        subplot(2,1,2); hold on; axis tight;

        means1_shuff = mean(cat(1, dec_data_means_shuff{:}),'omitnan');
        sem1_shuff = std(cat(1, dec_data_means_shuff{:}),'omitnan')/sqrt(size(cat(1, dec_data_means_shuff{:}),1)-1);
        shadedErrorBar_YS(params.(plot_var), means1_shuff, sem1_shuff, [0, 0, 0]);

        for n_reg = 1:numel(ops.regions_to_analyze)
            means1 = mean(dec_data_means{n_reg},'omitnan');
            sem1 = std(dec_data_means{n_reg},'omitnan')/sqrt(size(dec_data_means{n_reg},1)-1);
            shadedErrorBar_YS(params.(plot_var), means1, sem1, ops.cond_colors{n_reg});
        end
        legend([pl{:}, pls], [ops.regions_to_analyze; {'Shuff'}]);
        ylabel('performance')
        xlabel(strrep(plot_var,'_',' '))
        title(strrep(title_tag2,'_',' '))
        subplot(s1)
    else
        f_plot_dset_deets(dec_data_means, ops);
    end
    
end

end