function f_dv_plot_wf_tuning(app)

im_source = [2 7 9];

norm_within = 0;

ms_tag = app.ddata.mouse_tag{1};
reg_data = load(app.regdatapathEditField.Value);
ms_idx = strcmpi({reg_data.data_all.mouse_tag}, ms_tag);

map1 = reg_data.data_all(ms_idx).wf_mapping_im;
titl1 = reg_data.data_all(ms_idx).wf_mapping_title;

num_im = numel(map1);

[d1, d2] = size(map1{1});

figure; 
for n_im = 1:num_im
    subplot(2,5,n_im)
    im1 = imagesc(rot90(map1{n_im},1));
    axis equal tight
    im1.Parent.XTick = [];
    im1.Parent.YTick = [];
    title(titl1{n_im})
end
sgtitle(sprintf('%s; norm within = %d', ms_tag, norm_within), 'interpreter', 'none')

im_all = zeros(d1, d2, 3);
for n_sl = 1:3
    im_all(:,:,n_sl) = map1{im_source(n_sl)};
end 

if norm_within
    im_n = im_all;
    for n_chan = 1:3
        chan_slice = im_all(:,:,n_chan);
        im_n(:,:,n_chan) = im_all(:,:,n_chan)./max(chan_slice(:));
    end
else
    im_n = im_all./max(im_all(:));
end

figure;
im1 = imagesc(rot90(im_n*1.2,1));
axis equal tight
im1.Parent.XTick = [];
im1.Parent.YTick = [];
title(sprintf('%s; R-%s; G-%s; B-%s', ms_tag, titl1{im_source(1)}, titl1{im_source(2)}, titl1{im_source(3)}), 'interpreter', 'none')

end