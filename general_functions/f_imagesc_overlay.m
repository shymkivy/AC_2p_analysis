function f_imagesc_overlay(im1, im2)

imagecomb = zeros([size(im1),3]);
imagecomb(:,:,1) = im1/max(im1(:));
imagecomb(:,:,3) = im1/max(im1(:));
imagecomb(:,:,2) = im2/max(im2(:));
figure; imagesc(imagecomb)
title('fov registration check');

end