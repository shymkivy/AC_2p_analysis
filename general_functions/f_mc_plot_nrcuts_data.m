function f_mc_plot_nrcuts_data(cuts_data, title_tag)

% plot some nonrig stiff

    
num_planes = numel(cuts_data);

for n_pl = 1:num_planes
        
    if isfield(cuts_data{n_pl}, 'nr_corr_data')
        
        block_smooth = cuts_data{n_pl}.nr_corr_data.params.nonrigid_block_smooth;
        block_size = cuts_data{n_pl}.nr_corr_data.params.nonrigid_block_size;
        block_overlap = cuts_data{n_pl}.nr_corr_data.params.nonrigid_block_overlap;

        dsall_m = cuts_data{n_pl}.nr_corr_data.dsall_m;
        dsall_m_sm = f_smooth_movie_double(dsall_m, block_smooth);
        siz_nr = size(dsall_m);
        nr_2d = reshape(dsall_m, siz_nr(1)*siz_nr(2),siz_nr(3));
        nr_2d_sm = reshape(dsall_m_sm, siz_nr(1)*siz_nr(2),siz_nr(3));

%             dsall_mn = dsall_m_lr./std(dsall_m');
%             [U,S,V] = svd(nr_2d, 'econ');
%             n_comp = 10;
%             dsall_m_lr = U(:,1:n_comp)*S(1:n_comp,1:n_comp)*V(:,1:n_comp)';
%             
        figure; hold on;
        for n_rn = 1:min(siz_nr(1)*siz_nr(2),20)
            plot(nr_2d(n_rn,:) - 5*(n_rn-1))
            plot(nr_2d_sm(n_rn,:) - 5*(n_rn-1), 'k')
            %plot(dsall_m_lr(n_rn,:) - 5*(n_rn-1), 'g')
        end
        axis tight;
        ylabel(sprintf('blocks size = %d; overlap = %d', block_size, block_overlap))
        title(sprintf('%s, pl%d\nnonrigid raw vs block smooth [%.1f %.1f %.1f]',title_tag, n_pl, block_smooth(1), block_smooth(2), block_smooth(3)), 'interpreter', 'none');
    end
end

end