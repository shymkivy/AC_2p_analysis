function [thresh_coeffs, thresh_scores] = f_ens_get_thresh(firing_rate_ensemb, coeffs, scores, num_ens, params)
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf');
ensamble_extraction_thesh = f_get_param(params, 'ensamble_extraction_thesh', 'shuff'); % 'signal_z' 'shuff'
plot_stuff = f_get_param(params, 'plot_stuff', 0);

thresh_percent = 95;
shuff_rep = 50;
z_thresh = 2;

two_sided = logical(sum(coeffs(:)<0));
thresh_coeffs = zeros(num_ens,two_sided+1);
thresh_scores = zeros(num_ens,two_sided+1);

thresh_coeffs_z = zeros(num_ens,two_sided+1);
thresh_scores_z = zeros(num_ens,two_sided+1);


for n_comp = 1:num_ens
    factors = coeffs(:,n_comp);
    center1 = median(factors);
    z_fac = sqrt(sum((factors-center1).^2)/(numel(factors)-1));
    thresh_coeffs_z(n_comp,1) = center1 + z_thresh*z_fac;
    if two_sided
        thresh_coeffs_z(n_comp,2) = center1 - z_thresh*z_fac;
    end
end

for n_comp = 1:num_ens
    factors = scores(n_comp,:);
    center1 = median(factors);
    z_fac = sqrt(sum((factors-center1).^2)/(numel(factors)-1));
    thresh_scores_z(n_comp) = center1 + z_thresh*z_fac;
    if two_sided
        thresh_scores_z(n_comp,2) = center1 - z_thresh*z_fac;
    end
end

if strcmpi(ensamble_extraction_thesh, 'shuff')
    coeffs_shuff_all = cell(shuff_rep,1);
    scores_shuff_all = cell(shuff_rep,1);
    for n_rep = 1:shuff_rep
        firing_rate_ensemb_shuff = f_shuffle_data(firing_rate_ensemb, 'circ_shift');
        train_done = 0;
        while ~train_done
            try
                [dred_factors_shuff, ~] = f_dred_train2(firing_rate_ensemb_shuff, num_ens, ensamble_method, 0);
                train_done = 1;
            catch
                disp('Error train, will repeat');
            end
        end
        [coeffs_shuff, scores_shuff] = f_dred_get_coeffs(dred_factors_shuff);
        coeffs_shuff_all{n_rep} = coeffs_shuff;
        scores_shuff_all{n_rep} = scores_shuff';
    end
    coeffs_shuff_all1 = cat(1,coeffs_shuff_all{:});
    scores_shuff_all1 = cat(1,scores_shuff_all{:});
    if ~two_sided
        for n_comp = 1:num_ens
            thresh_coeffs(n_comp,1) = prctile(coeffs_shuff_all1(:,n_comp), thresh_percent);
            thresh_scores(n_comp,1) = prctile(scores_shuff_all1(:,n_comp), thresh_percent); 
        end
    else
        for n_comp = 1:num_ens
            thresh_coeffs(n_comp,:) = prctile(coeffs_shuff_all1(:,n_comp), [(100-thresh_percent)/2 thresh_percent+(100-thresh_percent)/2]);
            thresh_scores(n_comp,:) = prctile(scores_shuff_all1(:,n_comp), [(100-thresh_percent)/2 thresh_percent+(100-thresh_percent)/2]);
        end
    end
else
    thresh_coeffs = thresh_coeffs_z;
    thresh_scores = thresh_scores_z;
end

if plot_stuff
    figure;
    for n_comp = 1:num_ens
        subplot(num_ens,2,n_comp); hold on;
        stem(coeffs(:,n_comp))
        plot(ones(numel(coeffs(:,n_comp)),1)*thresh_coeffs_z(n_comp,:), '--r');
        if strcmpi(ensamble_extraction_thesh, 'shuff')
            plot(ones(numel(coeffs(:,n_comp)),1)*thresh_coeffs(n_comp,:), '--g');
        end
        axis tight;
        title(sprintf('Cells, comp %d', n_comp));
        xlabel('Cells')
    end
    for n_comp = 1:num_ens
        subplot(num_ens,2,n_comp+num_ens); hold on;
        stem(scores(n_comp,:))
        plot(ones(numel(scores(n_comp,:)),1)*thresh_scores_z(n_comp,:), '--r');
        if strcmpi(ensamble_extraction_thesh, 'shuff')
            plot(ones(numel(scores(n_comp,:)),1)*thresh_scores(n_comp,:), '--g');
        end
        axis tight;
    end
    if strcmpi(ensamble_extraction_thesh, 'shuff')
        legend('data', 'z thresh', 'shuff thresh');
    else
        legend('data', 'z thresh');
    end
end

end