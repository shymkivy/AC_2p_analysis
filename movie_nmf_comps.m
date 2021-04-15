
A = reshape(full(est.A),est.dims(1),est.dims(2),[]);

A3d = reshape(A, est.dims(1),est.dims(2),[]);

figure; imagesc(sum(A3d,3))

mov = reshape(full(est.A)*(est.C+est.YrA) + est.b'*est.f',est.dims(1),est.dims(2),[]);

bkg_im =mean(mov,3);
bkg_im = bkg_im/max(bkg_im(:));

figure; ax1 = imagesc(bkg_im); axis equal tight;

cells = [64, 2, 3, 31, 209, 75, 71, 93, 46, 21]; %rest 3
% cells = [2, 3, 4, 37, 46, 61]; % , 147, 182, 187, 104, 49%rest 4

rgbImage = ind2rgb(bkg_im, 'bone');

RGB = cat(3, bkg_im, bkg_im, bkg_im);

%figure; imagesc(RGB); axis equal tight;


cute_cells = A(:,:,cells);
cute_cells = cute_cells/max(cute_cells(:));

cute_cells2 = cute_cells;
for n_cell = 1:size(cute_cells,3)
    cute_cells2(:,:,n_cell) = imgaussfilt(cute_cells(:,:,n_cell), 1);
end

cute_cells_flat = reshape(cute_cells2, est.dims(1)*est.dims(2), []);

RGB_flat = reshape(RGB, [],3);
colormap1 = [.8 0 0];
im_loc = zeros(size(cute_cells_flat,1),3);
for n_cell = 1:numel(cells)
    temp = cute_cells_flat(:,n_cell)>0;
    for cc = 1:3
        im_loc(temp,cc) = cute_cells_flat(temp,n_cell)*colormap1(cc);
    end
end
sig_factor = 0.3/mean(im_loc(im_loc(:)>0));

imloc_idx = mean(im_loc,2);
imloc_idx = imloc_idx>0;

%RGB_flat(imloc_idx,:) = im_loc(imloc_idx,:);

RGB_flat = RGB_flat + im_loc;

RGB_flat3d = reshape(RGB_flat,est.dims(1),est.dims(2),3);

im_out = permute(mean(cute_cells,3), [2 1 3]);
im_rgb = permute(RGB_flat3d*sig_factor, [2 1 3]);

figure; imagesc(im_out)
imagesc(im_rgb); axis equal tight; axis off;
