function K = cosineKernel(vec1,vec2)

K = 1 - pdist2(vec1,vec2,'cosine');

end