function D = f_pdist2_YS(X, Y, distance_metric)

if strcmpi(distance_metric, 'hammilarity')
    [~, SI_hamm] = similarity_index(X,Y);
    D = (1 - SI_hamm);
else
    D = pdist2(X, Y, distance_metric);
end

end