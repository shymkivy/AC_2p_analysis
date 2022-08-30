% check moco resulst
clear
close all;

addpath([pwd '\s1_functions']);
addpath([pwd '\general_functions'])


average_dsbase =0;

%%
%ops.file_dir = 'F:\AC_data\caiman_data_dream3\preprocessing';
%ops.file_dir = 'F:\AC_data\caiman_data_echo\preprocessing';

params.dset_table_fpath = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\AC_data_list_all.xlsx';
%params.data_dir = 'F:\AC_data\caiman_data_dream';
%params.data_dir = 'D:\data\caiman_data_dream';
params.data_dir = 'F:\AC_data\caiman_data_missmatch';

% these have to be same as column names in excel
params.limit.dset_name =        '';
params.limit.experiment =       'missmatch';
params.limit.mouse_id =         'M10';
params.limit.mouse_tag =        '';
params.limit.dset_name =        '';
params.limit.FOV_num =          1;


%%
AC_data = f_s0_parse_tab_data(params);

%%
mouse_id_all = unique(AC_data.mouse_id, 'stable');

for n_ms1 = 1:numel(mouse_id_all)
    AC_data2 = AC_data(strcmpi(AC_data.mouse_id, mouse_id_all{n_ms1}),:);
    mouse_tag_all = unique(AC_data2.mouse_tag, 'stable');
    for n_ms2 = 1:numel(mouse_tag_all)
        AC_data3 = AC_data2(strcmpi(AC_data2.mouse_tag, mouse_tag_all{n_ms2}),:);
        
        mouse_fov_all = unique(AC_data3.FOV_num, 'stable');
        
        for n_fov = 1:numel(mouse_fov_all)
            AC_data4 = AC_data3(AC_data3.FOV_num == mouse_fov_all(n_fov),:);
            data_dir = [params.data_dir '\preprocessing\'];

            num_dsets = size(AC_data4,1);
            
            load_cuts_all = cell(num_dsets,1);
            for n_dset = 1:num_dsets
                flist = dir([data_dir, AC_data4.mouse_id{n_dset} '*' AC_data4.dset_name{n_dset} '*' AC_data4.mouse_tag{n_dset} '*h5cuts*']);
                load_cuts_all{n_dset} = load([data_dir flist.name]);
            end
            
            num_planes = numel(load_cuts_all{1}.cuts_data);
            leg_plane = cell(num_planes,1);
            for n_pl = 1:num_planes
                leg_plane{n_pl} = sprintf('plane %d', n_pl);
            end
            
            for n_dset = 1:num_dsets
                dset_fname = sprintf('%s_im%d_%s', AC_data4.mouse_id{n_dset}, AC_data4.im_num(n_dset), AC_data4.dset_name{n_dset});
                figure;
                for n_pl = 1:num_planes
                    subplot(num_planes, 1, n_pl); hold on;
                    ave_tr = load_cuts_all{n_dset}.cuts_data{n_pl}.ave_trace;
                    ave_tr = ave_tr - min(ave_tr);
                    ave_tr = ave_tr/max(ave_tr);
                    plot(ave_tr);
                    plot(load_cuts_all{n_dset}.cuts_data{n_pl}.vid_cuts_trace);
                end
                sgtitle(sprintf('vid cuts; %s', dset_fname), 'interpreter', 'none');
            end
            
            for n_pl = 1:num_planes
                figure;
                for n_dset = 1:num_dsets
                    dset_fname = sprintf('%s_im%d_%s', AC_data4.mouse_id{n_dset}, AC_data4.im_num(n_dset), AC_data4.dset_name{n_dset});
                    subplot(1, num_dsets, n_dset)
                    imagesc(load_cuts_all{n_dset}.cuts_data{n_pl}.image_target); axis equal tight;
                    title(sprintf('dset %s', dset_fname), 'interpreter', 'none')
                end
                sgtitle(sprintf('plane %d', n_pl))
            end
            
            moco_to_dset1 = unique(AC_data4.mc_to_dset, 'stable');
            idx1 = AC_data4.im_num == moco_to_dset1;
            ds_base_all = zeros(num_planes, num_dsets,2);
            ds_base_load = zeros(num_planes, num_dsets,2);
            corr_all = zeros(num_planes, num_dsets);
            im_reg_pre = cell(num_planes, num_dsets);
            % compute
            for n_pl = 1:num_planes
                im_target = load_cuts_all{idx1}.cuts_data{n_pl}.image_target;
                for n_dset = 1:num_dsets
                    im_target_dset = load_cuts_all{n_dset}.cuts_data{n_pl}.image_target;
                    [ds_base1, ~, ops_ds] = f_suite2p_reg_compute(im_target_dset, im_target);
                    ds_base_all(n_pl, n_dset, :) = ds_base1;
                    if isfield(load_cuts_all{n_dset}.cuts_data{n_pl}, 'ds_base')
                        ds_base_load(n_pl, n_dset, :) = load_cuts_all{n_dset}.cuts_data{n_pl}.ds_base;
                    else
                        ds_base_load(n_pl, n_dset, :) = [0, 0];
                    end
                    corr_all(n_pl, n_dset) = ops_ds{1}.CorrFrame;
                    im_reg_pre{n_pl, n_dset} = im_target_dset;
                end
            end
            % apply
            im_reg_post = cell(num_planes, num_dsets);
            for n_pl = 1:num_planes
                im_target = load_cuts_all{idx1}.cuts_data{n_pl}.image_target;
                for n_dset = 1:num_dsets
                    im_target_dset = load_cuts_all{n_dset}.cuts_data{n_pl}.image_target;
                    if average_dsbase
                        ds_base1 = squeeze(mean(ds_base_all(:, n_dset, :),1))';
                    else
                        ds_base1 = squeeze(ds_base_all(n_pl, n_dset, :))';
                    end
                    im_reg_post{n_pl, n_dset} = f_suite2p_reg_apply(im_target_dset, ds_base1);
                end
            end
            
            
            SI_pre = cell(num_planes,1);
            SI_post = cell(num_planes,1);
            for n_pl = 1:num_planes
                im_reg2 = cat(3,im_reg_pre{n_pl,:});
                [d1, d2, ~] = size(im_reg2);
                im_reg3 = reshape(im_reg2, d1*d2, num_dsets);
                SI_pre{n_pl} = 1 - pdist2(im_reg3', im_reg3', 'correlation');
                im_reg2 = cat(3,im_reg_post{n_pl,:});
                im_reg3 = reshape(im_reg2, d1*d2, num_dsets);
                SI_post{n_pl} = 1 - pdist2(im_reg3', im_reg3', 'correlation');
            end
            if 0
                si_cat_post = cat(1, SI_post{:});
                clim_si_post = [min(si_cat_post(:)) max(si_cat_post(:))];
                si_cat_diff = cat(1, SI_post{:}) - cat(1, SI_pre{:});
                clim_si_diff = [min(si_cat_diff(:)) max(si_cat_diff(:))];
                for n_pl = 1:num_planes
                    figure; imagesc(SI_post{n_pl});
                    title(sprintf('dset SI post plane %d', n_pl))
                    caxis(clim_si_post);

                    figure; imagesc(SI_post{n_pl} - SI_pre{n_pl});
                    title(sprintf('dset SI post - pre plane %d', n_pl))
                    caxis(clim_si_diff);
                end
            end
            
            if 0
                for n_pl = 1:num_planes
                    figure;
                    for n_dset = 1:num_dsets
                        dset_fname = sprintf('%s_im%d_%s', AC_data4.mouse_id{n_dset}, AC_data4.im_num(n_dset), AC_data4.dset_name{n_dset});
                        subplot(1, num_dsets, n_dset)
                        imagesc(im_reg_post{n_pl,n_dset}); axis equal tight;
                        title(sprintf('dset %s', dset_fname), 'interpreter', 'none')
                    end
                    sgtitle(sprintf('post reg; plane %d', n_pl))
                end
            end
            
            % plot reg in same plot
            num_plots = ceil((num_dsets-1)/2);
            for n_pl = 1:num_planes
                for n_plt = 1:num_plots
                    num_col = min(num_dsets - 1 - 2*(n_plt-1),2);
                    im3_pre = zeros(d1,d2,3);
                    im3_pre(:,:,1) = im_reg_pre{n_pl,1};
                    im3_post = zeros(d1,d2,3);
                    im3_post(:,:,1) = im_reg_post{n_pl,1};
                    for n_col = 1:num_col
                        im3_pre(:,:,n_col+1) = im_reg_pre{n_pl,(n_plt-1)*2 + n_col+1};
                        im3_post(:,:,n_col+1) = im_reg_post{n_pl,(n_plt-1)*2 + n_col+1};
                    end
                    im3_pre = im3_pre/max(im3_pre(:));
                    im3_post = im3_post/max(im3_post(:));
                    
                    figure;
                    subplot(1,2,1);
                    imagesc(im3_pre*2); title('pre')
                    subplot(1,2,2);
                    imagesc(im3_post*2); title('post')
                    sgtitle(sprintf('plane %d, dsets 1, %d, %d', n_pl, (n_plt-1)*2 + 1 +1, (n_plt-1)*2 + 2 +1))
                end
            end
                    
            figure; 
            subplot(2,1,1);
            plot(ds_base_all(:,:,1), '-o');
            ylabel('y offset');
            subplot(2,1,2);
            plot(ds_base_all(:,:,2), '-o')
            ylabel('x offset');
            xlabel('num plane');
            sgtitle('dsbase across planes')
            
            if 1
                for n_dset = 1:num_dsets
                    
                    dset_fname = sprintf('%s_im%d_%s; pl%d', AC_data4.mouse_id{n_dset}, AC_data4.im_num(n_dset), AC_data4.dset_name{n_dset}, n_pl);   
                    cuts_data1 = load_cuts_all{n_dset}.cuts_data;
                    f_mc_plot_cuts_data(cuts_data1, dset_fname);
                    
                    f_mc_plot_nrcuts_data(cuts_data1, dset_fname)
                    
                end
            end
            
            if 1
                data_dir = [params.data_dir '\movies\'];
                ave_im_all = cell(num_planes, num_dsets);
                for n_pl = 1:num_planes
                    for n_dset = 1:num_dsets
                        if num_planes > 1
                            fsearch = sprintf('%s*%s*%s*%s*_pl%d.h5', data_dir, AC_data4.mouse_id{n_dset},  AC_data4.dset_name{n_dset}, AC_data4.mouse_tag{n_dset}, n_pl);
                        else
                            fsearch = sprintf('%s*%s*%s*%s.h5', data_dir, AC_data4.mouse_id{n_dset},  AC_data4.dset_name{n_dset}, AC_data4.mouse_tag{n_dset});
                        end
                        flist = dir(fsearch);
                        temp_fpath = [data_dir '\' flist.name];
                        if numel(flist)
                            Y_temp = single(h5read(temp_fpath,'/mov'));
                            ave_im_all{n_pl, n_dset} = mean(Y_temp,3);
                        end
                    end
                    for n_plt = 1:num_plots
                        if ~isempty(ave_im_all{n_pl, n_dset})
                            num_col = min(num_dsets - 1 - 2*(n_plt-1),2);
                            im3_pre = zeros(d1,d2,3);
                            im3_pre(:,:,1) = ave_im_all{n_pl,1};
                            for n_col = 1:num_col
                                im3_pre(:,:,n_col+1) = ave_im_all{n_pl,(n_plt-1)*2 + n_col+1};
                            end
                            im3_pre = im3_pre/max(im3_pre(:));

                            figure;
                            imagesc(im3_pre*2); 
                            if num_col > 1
                                title(sprintf('from h5; plane %d, dsets 1, %d, %d', n_pl, (n_plt-1)*2+2, (n_plt-1)*2+3));
                            else
                                title(sprintf('from h5; plane %d, dsets 1, %d', n_pl, (n_plt-1)*2+2));
                            end
                        end
                    end
                    
                end
            end
            
        end
    end
end
disp('Done')
