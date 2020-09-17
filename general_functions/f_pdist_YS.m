function D = f_pdist_YS(X, distance_metric)

if strcmpi(distance_metric, 'hammilarity')
    [~, SI_hamm] = similarity_index(X,X);
    min_d = min(size(SI_hamm));
    zero_diag_mat = 1-diag(ones(min_d,1));
    D = (1 - SI_hamm).*zero_diag_mat;
else
    D = squareform(pdist(X, distance_metric));
end

end