num_d = 10;
num_pts = 1000;
x = zeros(num_d, num_pts);

for n_x = 1:num_d
    x(n_x,:) = randn(1000,1);
end

mean_vec = mean(x,2);

dist1 = pdist2(mean_vec', x', 'cosine'); % euclidean cosine

figure; histogram(dist1)