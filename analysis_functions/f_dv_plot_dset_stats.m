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
    mmn_all2 = mmn_all(:);
    idx_mid = and(mmn_all2 > 1, mmn_all2<9);
    
    fprintf('%s; stim diff mean %.2f %s, std %.2f %s\n', data.paradigm{1}, mean(stim_val_diff), units1, std(stim_val_diff), units1);
    fprintf('min MMN stim %d; max MMN stim %d\n', min(mmn_all(:)), max(mmn_all(:)));

    figure();
    histogram(mmn_all2, 'BinEdges', (0:10)+0.5);
    title(sprintf('oddball stimuli used; %s', data(1,:).paradigm{1}), 'interpreter', 'none');
    fprintf('%.2f%% including 2-8\n', sum(idx_mid)/numel(idx_mid)*100);
end


rec_day = data.im_recording_day/7;
age = data.Cran_age;

num_tun_uq = zeros(num_dsets,1);
num_tun_all = zeros(num_dsets,10);
for n_dset = 1:num_dsets
    num_tun_uq(n_dset) = sum(logical(sum(data(n_dset,:).stats{1}.peak_resp_cells(:,1:10),2)),1);
    num_tun_all(n_dset,:) = sum(data(n_dset,:).stats{1}.peak_resp_cells(:,1:10));
end

use_data = and(and(~isnan(age), rec_day<2.5), num_tun_uq<100);

age2 = age(use_data) + rec_day(use_data);
num_tun_uq2 = num_tun_uq(use_data);
num_tun_all2 = num_tun_all(use_data,:);


mtag = data.mouse_tag(use_data);
unique(mtag)


x_fit = 13:0.01:18;
fitob = fit(age2,num_tun_uq2,'poly1');
y_fit = fitob(x_fit);

figure(); hold on;
plot(age2, num_tun_uq2, '.');
plot(x_fit, y_fit, color='k');
ylabel('number of tuned cells');
xlabel('age at experiment (weeks)');
title('age vs tuning counts');
axis padded

col = jet(10);

figure(); hold on;
for n_fr = 1:10
    fitob = fit(age2,num_tun_all2(:,n_fr),'poly1');
    y_fit = fitob(x_fit);
    plot(age2, num_tun_all2(:,n_fr), '.', color=col(n_fr,:));
    plot(x_fit, y_fit, color=col(n_fr,:));
end
ylabel('number of tuned cells');
xlabel('age at experiment (weeks)');
title('age vs tuning counts');
axis padded

end 
