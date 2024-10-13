function f_dv_plot_dset_stats(app)

data = app.data;

num_cells2 = cell2mat(data.num_cells_pl);
num_cells3 = sum(num_cells2,2);


num_dsets = numel(num_cells3);
num_cells_mean = mean(num_cells3);
num_cells_std = std(num_cells3);
num_cells_sem = num_cells_std/sqrt(numel(num_cells3)-1);

fprintf('%d datasets; mean %.2f std %.2f sem %.2f cells per dset\n', num_dsets, num_cells_mean, num_cells_std, num_cells_sem);

if sum(strcmpi(data(1,:).paradigm, {'tone_mmn', 'fg_mmn'}))
    mmn_all = cat(1,data.MMN_freq{:});
    num_dsets = size(data,1);
    
    stim_num_diff = zeros(num_dsets,1);
    num_stim = zeros(num_dsets,1);
    stim_diff = zeros(num_dsets,1);
    if strcmpi(data.paradigm{1}, 'fg_mmn') % circular stim space
    
        
        ori_all = cell(num_dsets,1);
        for n_dset = 1:num_dsets
            data2 = data(n_dset,:).proc_data{1}.stim_params;
            if isfield(data2, 'grating_angles')
                num_stim(n_dset) = numel(data2.grating_angles);
                ori_all{n_dset} = data2.grating_angles;
            elseif isfield(data2, 'ops')
                num_stim(n_dset) = numel(data2.ops.grating_angles);
                ori_all{n_dset} = data2.ops.grating_angles;
            end
            stim_diff(n_dset) = abs(mean(diff(ori_all{n_dset})));
        end
    
        for n_dset = 1:num_dsets
            [~, mx_ind] = max(mmn_all(n_dset,:));
    
            diff1 = mmn_all(n_dset, mx_ind) - mmn_all(n_dset, 3-mx_ind);
            diff2 = mmn_all(n_dset, 3-mx_ind)+num_stim(n_dset) - mmn_all(n_dset, mx_ind);
    
            stim_num_diff(n_dset) = min([diff1, diff2]);
        end
        
        stim_val_diff = stim_num_diff.*stim_diff/pi;
    
        units1 = 'pi rad';
        
        x1 = real(exp(1i*ori_all{1}*2));
        y1 = imag(exp(1i*ori_all{1}*2));
        
        scale1 = 20;
        figure(); hold on;
        for n_dset = 1:num_dsets
            plot(x1(mmn_all(n_dset,:)) + (rand(1,2)-0.5)/scale1, y1(mmn_all(n_dset,:)) + (rand(1,2)-0.5)/scale1, 'o-k');
        end
        plot(x1, y1, 'bo', linewidth=2);
        axis equal padded;
    
    elseif strcmpi(data.paradigm{1}, 'tone_mmn') % linear stim space
        for n_dset = 1:num_dsets
            data2 = data(n_dset,:).proc_data{1}.stim_params;
            num_stim(n_dset) = data2.num_freqs;
            stim_diff(n_dset) = data2.increase_factor;
        end
    
        for n_dset = 1:num_dsets
            [~, mx_ind] = max(mmn_all(n_dset,:));
    
            stim_num_diff(n_dset) = mmn_all(n_dset, mx_ind) - mmn_all(n_dset, 3-mx_ind);
        end
        
        stim_val_diff = stim_num_diff.*stim_diff;
        
        units1 = 'octaves';
    end
    
    fprintf('%s; stim diff mean %.2f %s, std %.2f %s\n', data.paradigm{1}, mean(stim_val_diff), units1, std(stim_val_diff), units1);
end


end 
