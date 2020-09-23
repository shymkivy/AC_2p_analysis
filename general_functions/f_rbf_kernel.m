function dist_out = f_rbf_kernel(vec1, vec2, gamma)

if ~exist('gamma', 'var') || isempty(gamma)
    gamma = 2;
end

%d1 = size(vec1,1);
%d2 = size(vec2,1);
%mat1 = repmat(reshape(vec1,d1, 1, []),1,d2,1);
%mat2 = repmat(reshape(vec2,1, d2, []),d1,1,1);
%dist_out = exp(-sum((reshape(vec1,d1, 1, []) - reshape(vec2,1, d2, [])).^2,3));

dist_out = exp(-pdist2(vec1,vec2, 'euclidean').^2./(2*gamma^2));

%figure; imagesc(dist_out)

end