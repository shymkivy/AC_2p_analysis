close all
clear

%%
addpath([pwd '\s1_functions']);
addpath([pwd '\general_functions'])

%%
%ops.file_dir = 'F:\AC_data\caiman_data_dream3\preprocessing';
%ops.file_dir = 'F:\AC_data\caiman_data_echo\preprocessing';

params.dset_table_fpath = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\AC_data_list_all.xlsx';
%params.data_dir = 'F:\AC_data\caiman_data_dream3';
params.data_dir = 'D:\data\caiman_data_dream';

% these have to be same as column names in excel
params.limit.dset_name =        '';
params.limit.experiment =       'dream';
params.limit.mouse_id =         'M125';
params.limit.mouse_tag =        '4_26_22';
params.limit.dset_name =        '';

roi_half_size = 6;

%%
AC_data = f_s0_parse_tab_data(params);

%%
mouse_id_all = unique(AC_data.mouse_id, 'stable');

% work plane by plane

num_dsets = numel(AC_data.mouse_id);
num_planes = max(AC_data.mpl);

ave_frames = cell(num_planes, num_dsets);
std_frames = cell(num_planes, num_dsets);
num_frames = zeros(num_planes, num_dsets);

proc_data = cell(num_dsets,1);
for n_dset = 1:num_dsets
    flist = dir([params.data_dir '\preprocessing\*' AC_data.dset_name{n_dset} '*' AC_data.mouse_tag{n_dset} '*processed_data*']);
    if numel(flist.name)
        proc_data{n_dset} = load([params.data_dir '\preprocessing\' flist.name]);
    end
end


cell_mn2 = cell(num_planes,1);
cell_spat2 = cell(num_planes,1);
cell_ca2 = cell(num_planes,1);
cell_run_id2 = cell(num_planes,1);
cell_pl_id = {};

for n_pl = 1:num_planes
    Y = cell(num_dsets,1);
    for n_dset = 1:num_dsets
        flist = dir([params.data_dir '\movies\*' AC_data.dset_name{n_dset} '*' AC_data.mouse_tag{n_dset} '*_pl' num2str(n_pl) '*']);
        if numel(flist.name)
            Y_temp = h5read([params.data_dir '\movies\' flist.name], '/mov');
            [d1, d2, ~] = size(Y_temp);
            vid_cuts_trace = logical(proc_data{n_dset}.data.file_cuts_params{n_pl}.vid_cuts_trace);
            Y{n_dset} = zeros(d1, d2, numel(vid_cuts_trace), 'uint16');
            Y{n_dset}(:,:,vid_cuts_trace) = Y_temp;
            
            ave_frames{n_pl, n_dset} = mean(single(Y{n_dset}),3);
            std_frames{n_pl, n_dset} = std(single(Y{n_dset}),0,3);
            num_frames(n_pl, n_dset) = size(Y{n_dset},3);
        end
        
    end
    Y_all = cat(3, Y{:});
    clear Y Y_temp; 
    
    im1 = cat(2,std_frames{n_pl, :});
    
    figure;
    imagesc(im1);
    axis equal tight;
    title(sprintf('%s; plane %d; dset 1 - %d', AC_data.mouse_id{n_dset}, n_pl, num_dsets), 'interpreter', 'none');
    
    
    ave_frame_all = mean(single(Y_all),3);
    std_frame_all = std(single(Y_all), 0, 3);
    
    Y_alln = f_s4_remove_mov_bkg(Y_all, 1);
 
    title_all = {'target_cell', 'coactivated_cell', 'other'};
    color_all = {'r', 'g', 'b'};
    stim_coact_other_cells_mn = cell(0,3);
    
    cell_mn = {};
    cell_spat = {};
    cell_ca = {};
    cell_run_id = {};
    
    f2 = figure;
    f1 = figure; hold on;
    im1 = imagesc(std_frame_all); axis equal tight;
    im1.Parent.YDir = 'reverse';
    
    for n_run = 1:3
        n_try = 1;
        ans1 = 1;
        
        while ans1
            
            figure(f1);
            title(sprintf('plane %d; add %s cell? (y/n)', n_pl, title_all{n_run}), 'interpreter', 'none');
            ans1 = ask_yes_no_fig();
            if ans1
                n_try = n_try + 1;
                title('Click on center of roi', 'interpreter', 'none');
                [n, m] = ginput(1);
                mr = round(m);
                nr = round(n); 
                % idx1[m1, m2, n1, n2]
                idx1 = [(mr-roi_half_size), (mr+roi_half_size); (nr-roi_half_size), (nr+roi_half_size)];
                
                % get mean trace and roi
                roi1 = Y_all(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2),:);
                roi1n = Y_alln(idx1(1,1):idx1(1,2), idx1(2,1):idx1(2,2),:);
                trace1 = squeeze(mean(mean(roi1,1),2));
                trace1 = trace1 - mean(trace1);
                
                [dr1, dr2, Tr] = size(roi1n);
                roi1n2d = reshape(roi1n, dr1*dr2, Tr);
                
                % get svd version
                %U*S*V'
                [U,S,V] = svd(single(roi1n2d),'econ');
                %data_1d = U(:,1)*S(1,1)*V(:,1)';
                spat1 = reshape(U(:,1)*S(1,1)*mean(V(:,1)),dr1, dr2);
                spat_norm = norm(spat1(:));
                trace1n = mean(U(:,1))*S(1,1)*V(:,1);
                spat1 = spat1/spat_norm;
                trace1n = trace1n/spat_norm;
                
                if mean(mean(spat1(:,[1, round(2*roi_half_size+1)]))) < 0
                    spat1 = spat1 * -1;
                end
                
                %trace1nd = f_smooth_dfdt3(trace1n', 1, 2);
                figure(f2);
                subplot(2,2,1);
                imagesc(mean(roi1,3))
                title('mean roi frame')
                subplot(2,2,2);
                imagesc(spat1)
                title('svd extracted roi')
                subplot(2,2,3:4); hold on
                plot(trace1/norm(trace1));
                plot(trace1n/norm(trace1n));
                axis tight;
                title('mean roi trace')
                legend('mean', 'svd extracted')
                sgtitle(sprintf('%s, plane %d; try %d', AC_data.mouse_id{1}, n_pl, n_try))
                
                
                
                figure(f1);
                r1 = rectangle('Position', [idx1(2,1), idx1(1,1), abs(diff(idx1(1,:))), abs(diff(idx1(1,:)))], 'EdgeColor', color_all{n_run});
                title('is it good?', 'interpreter', 'none');
                ans2 = ask_yes_no_fig();
                
                if ans2
                    cell_mn = [cell_mn; mr nr];
                    cell_spat = [cell_spat; spat1];
                    cell_ca = [cell_ca, trace1n'];
                    cell_run_id = [cell_run_id; n_run];
                else
                    delete(r1)
                end
            end
        end
    end
    
    title('Done', 'interpreter', 'none');
    
    cell_mn2{n_pl} = cat(1, cell_mn{:});
    cell_spat2{n_pl}  = cat(3,cell_spat{:});
    cell_ca2{n_pl}  = cat(1, cell_ca{:});
    cell_run_id2{n_pl}  = cat(1,cell_run_id{:});
    cell_pl_id = [cell_pl_id; ones(numel(cell_run_id),1)*n_pl];
end    

%%
data.AC_data = AC_data;
data.cell_mn = cell_mn2;
data.cell_spat = cell_spat2;
data.cell_ca = cell_ca2;
data.cell_run_id = cell_run_id2;
data.cell_pl_id = cat(1,cell_pl_id{:});

data.ave_frames = ave_frames;
data.std_frames = std_frames;
data.num_frames = num_frames;

data.proc_data = proc_data;

%%  
save_name = [params.data_dir '\preprocessing\' AC_data.mouse_id{n_dset} '_' AC_data.dset_name{n_dset} '_' AC_data.mouse_tag{n_dset} '_manual_roi.mat'];
    
save(save_name, 'data', 'params');
