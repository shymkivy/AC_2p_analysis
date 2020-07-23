function SI = similarity_index(vec1, vec2)
% calculates Similarity Index
% we use this instead cosine of angle because this accounts for magnitude
% and cosine doesn't
% SI = (Ca * Cb) / ((||Ca||^2 + ||Cb||^2)/2) - accounts for magnitude diff
% best for highly skewed distributions
% cos(ang) = (Ca * Cb) / (||Ca|| * ||Cb||) - ok for binary
    dot_prod = vec1*vec2';
    magnitudes1 = ones(size(vec1,1))*diag(diag(vec1*vec1'));
    magnitudes2 = ones(size(vec2,1))*diag(diag(vec2*vec2'));
    SI = dot_prod./((magnitudes1' + magnitudes2)/2);
end