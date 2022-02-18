function f_dv_register_cells(app)

fov_idx = [4,7];

bkg_all = cell(numel(fov_idx),1);
for n_idx = 1:numel(fov_idx)
    data = app.data(fov_idx(n_idx),:);
    est = data.OA_data{1}.est;
    
    bkg = reshape(mean(est.f)*est.b, [256 256]);
    figure; imagesc(bkg); title(num2str(fov_idx(n_idx)));
    bkg_all{n_idx} = bkg;
end

[d11, d21] = size(bkg_all{1});
[d12, d22] = size(bkg_all{2});

df1 = d11+d12-1;
df2 = d21+d22-1;

four1 = fft2(bkg_all{1}, df1, df2);
four2 = fft2(rot90(bkg_all{2},2), df1, df2); % flip in both directions

back1 = ifft2(four1.*four2);

[~, idx] = max(back1(:));

[m1, n1] = ind2sub([df1,df2], idx);

shiftm = m1 - d11;
shiftn = n1 - d12;

A_all = cell(numel(fov_idx),1);
A_all3 = cell(numel(fov_idx),1);
for n_idx = 1:numel(fov_idx)
    data = app.data(fov_idx(n_idx),:);
    est = data.OA_data{1}.est;
    proc = data.OA_data{1}.proc;
    A_full = reshape(full(est.A(:,proc.comp_accepted)), d11, d12, []);
    A_flat = sum(A_full,3);
    %figure; imagesc(A_flat);
    %title(num2str(fov_idx(n_idx)));
    A_all3{n_idx} = A_full;
    A_all{n_idx} = A_flat/max(A_flat(:))*2;
end

A_rgb = cat(3,A_all{:}, zeros(256, 256));

figure; imagesc(A_rgb)


A_all3_n = A_all3;
for n_dset = 1:numel(A_all3)
    A = A_all3{n_dset};
    A_n = A;
    for n_cell = 1:size(A_n,3)
        A_cell = A(:,:,n_cell);
        mask1 = A_cell>0;
        A_cell_n = A_cell;
        A_cell_n(mask1) = A_cell_n(mask1) - mean(A_cell_n(mask1));
        A_cell_n(mask1) = A_cell_n(mask1)/std(A_cell_n(mask1));
        A_n(:,:,n_cell) = A_cell_n;
    end
    A_all3_n{n_dset} = A_n;
end

dset1 = A_all3_n{1};
dset2 = A_all3_n{1};

num_cells1 = size(dset1,3);
num_cells2 = size(dset2,3);

c_vals = zeros(num_cells1, num_cells2);
c_dist = zeros(num_cells1, num_cells2);

for n_cell1 = 1:num_cells1
    for n_cell2 = 1:num_cells2
        cell1 = dset1(:,:,n_cell1);
        cell2 = dset2(:,:,n_cell2);
        convvv = ifft2(fft2(cell1, 511, 511).*fft2(rot90(cell2,2), 511, 511));
        
        [c_vals(n_cell1, n_cell2), max_idx] = max(convvv(:));
        [m1c, n1c] = ind2sub(size(convvv), max_idx);
        
        c_dist(n_cell1, n_cell2) = sqrt((m1c-d11)^2 + (n1c-d12)^2);
    end
end

dset1 = A_all3_n{1};
dset2 = A_all3_n{2};

num_cells1 = size(dset1,3);
num_cells2 = size(dset2,3);

c_vals2 = zeros(num_cells1, num_cells2);
c_dist2 = zeros(num_cells1, num_cells2);

for n_cell1 = 1:num_cells1
    for n_cell2 = 1:num_cells2
        cell1 = dset1(:,:,n_cell1);
        cell2 = dset2(:,:,n_cell2);
        convvv = ifft2(fft2(cell1, 511, 511).*fft2(rot90(cell2,2), 511, 511));
        
        [c_vals2(n_cell1, n_cell2), max_idx] = max(convvv(:));
        [m1c, n1c] = ind2sub(size(convvv), max_idx);
        
        c_dist2(n_cell1, n_cell2) = sqrt((m1c-d11)^2 + (n1c-d12)^2);
        
        if c_dist2(n_cell1, n_cell2) < 5
            c_rgb = cat(3,cell1/max(cell1(:))*1.2, circshift(cell2, 20, 2)/max(cell2(:))*1.2, zeros(256, 256));
            figure; imagesc(c_rgb);
            title(sprintf('cell %d, cell %d, c val %.2f', n_cell1, n_cell2, c_vals2(n_cell1, n_cell2)))
        end
    end
end



figure; imagesc(c_vals)
figure; imagesc(c_dist)

figure; hold on;
plot(c_vals(:), c_dist(:), 'o', 'color', [.4 .4 .4])
plot(c_vals2(:), c_dist2(:), 'o', 'color', 'b')


cell1 = A_all3_n{1}(:,:,1);
cell2 = A_all3_n{1}(:,:,2);
convvv = ifft2(fft2(cell1, 511, 511).*fft2(rot90(cell2,2), 511, 511));

[max_val, max_idx] = max(convvv(:));
[m1c, n1c] = ind2sub(size(convvv), max_idx);





m1 - m2
n1 - n2

m1*2-1


figure; imagesc(cell1)
figure; imagesc(convvv)

figure; imagesc(c2)

norm(A_cell(mask1)-mean(A_cell(mask1)), 'fro')

var(A_cell(mask1))

A_cell2 = A(:,:,1);
A_cell2_n = A_cell2/norm(A_cell2, 'fro');

max(max(ifft2(fft2(A_cell_n, 511, 511).*fft2(A_cell_n, 511, 511))))


figure; imagesc(ifft2(fft2(A_n(:,:,1), 511, 511).*fft2(A_n(:,:,3), 511, 511)))

end