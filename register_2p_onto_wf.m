

wf_im_path = 'G:\data\AC\mapping\5_17_20_mapping\testing\Average Image AC_mapping1_5_17_20.fig';
twop_path = 'G:\data\AC\2p\5_17_20_im\AAF_surface-001\AAF_surface-001_Cycle00001_Ch2_000001.ome.tif';

fig1 = open(wf_im_path);
wf_im = fig1.Children.Children.CData;
twop_im = imread(twop_path);

wf_im2 = wf_im - min(wf_im(:));
wf_im2 = wf_im2/max(wf_im2(:));
%wf_im2 = wf_im2*100;

twop_im2 = double(twop_im);
twop_im2 = twop_im2 - min(twop_im2(:));
twop_im2 = twop_im2/max(twop_im2(:));
%twop_im2 = twop_im2*256;

imwrite(wf_im2,'reg_fixed_im.png','PNG');
imwrite(twop_im2, 'reg_moving_im.png','PNG'); %,parula(1000)

[xy2p2, xywf2] = cpselect('reg_moving_im.png','reg_fixed_im.png','Wait',true);

tform = fitgeotrans(xy2p2,xywf2,'affine');

movingRegistered = imwarp(twop_im2,tform,'OutputView',imref2d(size(wf_im)));
%figure; imagesc(movingRegistered)

comb_im = wf_im2;
comb_im(movingRegistered>0) = movingRegistered(movingRegistered>0);
figure; 
subplot(1,2,1);
imagesc(wf_im2);axis equal tight;
subplot(1,2,2);
imagesc(comb_im);axis equal tight;


% 
% f_2p = figure;
% imagesc(twop_im); hold on;
% axis equal tight;
% 
% f_wf = figure; 
% imagesc(wf_im); hold on;
% axis equal tight;
% 
% num_points = 3;
% xy2p = zeros(num_points,2);
% xywf = zeros(num_points,2);
% 
% for n_pt = 1:num_points
%     figure(f_2p);
%     title(['Select point ' num2str(n_pt)]);
%     [x1, y1] = ginputColor(1);
%     title('Other plot');
%     xy2p(n_pt,:) = round([x1 y1]);
%     plot(xy2p(n_pt,1),xy2p(n_pt,2), 'om');
%     text(xy2p(n_pt,1)+2,xy2p(n_pt,2),num2str(n_pt))
%     figure(f_wf);
%     title(['Select corresponding point ' num2str(n_pt)]);
%     [x1, y1] = ginputColor(1);
%     title('Other plot');
%     xywf(n_pt,:) = round([x1 y1]);
%     plot(xywf(n_pt,1),xywf(n_pt,2), 'om');
%     text(xywf(n_pt,1)+2,xywf(n_pt,2),num2str(n_pt))
% end
