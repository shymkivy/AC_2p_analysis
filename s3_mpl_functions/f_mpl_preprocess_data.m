function data = f_mpl_preprocess_data(data, ops)

disp('Preprocessing...');

num_dsets = numel(data.area);
if ops.waitbar
    wb = f_waitbar_initialize([], 'Loading data...');
end
for n_dset = 1:num_dsets
    temp_proc_data = data.proc_data{n_dset,1};
    frame_period_sec = temp_proc_data.frame_data.volume_period/1000;
    
    data.trial_types{n_dset} = temp_proc_data.trial_types;
    if isfield(temp_proc_data.stim_params, 'MMN_freq')
        data.MMN_freq{n_dset} = temp_proc_data.stim_params.MMN_freq;
    else
        data.MMN_freq{n_dset} = temp_proc_data.MMN_orientations;
    end
    
    data.trial_window{n_dset} = struct();
    data.trial_window{n_dset}.trial_window_t = (ceil(ops.trial_window(1)/frame_period_sec):floor(ops.trial_window(2)/frame_period_sec))*frame_period_sec;
    %data.trial_window_t_long{n_dset} = (ceil(ops.trial_window_long(1)/frame_period_sec):floor(ops.trial_window_long(2)/frame_period_sec))*frame_period_sec;
    data.trial_window{n_dset}.trial_num_baseline_resp_frames = [sum(data.trial_window{n_dset}.trial_window_t<=0) sum(data.trial_window{n_dset}.trial_window_t>0)];     
    %data.trial_num_baseline_resp_frames_long{n_dset} = [sum(data.trial_window_t_long{n_dset}<=0) sum(data.trial_window_t_long{n_dset}>0)];  
    data.trial_window{n_dset}.baseline_window_frames = data.trial_window{n_dset}.trial_window_t<=0;
    data.trial_window{n_dset}.resp_window_frames = and(data.trial_window{n_dset}.trial_window_t>ops.resp_window_time(1),data.trial_window{n_dset}.trial_window_t<ops.resp_window_time(2));       
    data.trial_window{n_dset}.resp_window_t = data.trial_window{n_dset}.trial_window_t(data.trial_window{n_dset}.resp_window_frames);
    data.trial_window{n_dset}.onset_window_frames = and(data.trial_window{n_dset}.trial_window_t>ops.onset_window(1),data.trial_window{n_dset}.trial_window_t<ops.onset_window(2)); 
    data.trial_window{n_dset}.offset_window_frames = and(data.trial_window{n_dset}.trial_window_t>ops.offset_window(1),data.trial_window{n_dset}.trial_window_t<ops.offset_window(2)); 

    cell_plane_indx_pl = cell(data.num_planes(n_dset),1);
    
    % pull out data
    for n_pl = 1:data.num_planes(n_dset)
        temp_OA_data = data.OA_data{n_dset,n_pl};
        temp_proc_data = data.proc_data{n_dset};
        % extra cSNR threshold
        SNR_accept = temp_OA_data.proc.SNR2_vals >= ops.extra_SNR_thresh;
        accept_cell = and(SNR_accept,temp_OA_data.proc.comp_accepted);           
        traces_raw_cut = temp_OA_data.est.C(accept_cell,:)+temp_OA_data.est.YrA(accept_cell,:);

        try
            firing_rate_cut = if_get_deconvolved_data(temp_OA_data,accept_cell, ops,temp_proc_data.frame_data.volume_period_ave);
        catch
            error([data.experiment{n_dset} ' deconvolution incomplete'])
        end
        cuts_trace = temp_proc_data.file_cuts_params{n_pl}.vid_cuts_trace;

        if ops.normalize_firing_rate
            firing_rate_cut = (firing_rate_cut - min(firing_rate_cut,[],2))./max((firing_rate_cut - min(firing_rate_cut,[],2)),[],2);
        end

        if ops.signal_extra_smooth_sig
            % additional smooth 
            kernel_sigma_frames = ops.signal_extra_smooth_sig/frame_period_sec;
            kernel_half_size = ceil(sqrt(-log(0.05)*2*kernel_sigma_frames^2));
            gaus_win = -kernel_half_size:kernel_half_size;
            gaus_kernel = exp(-((gaus_win).^2)/(2*kernel_sigma_frames^2));
            gaus_kernel = gaus_kernel/sum(gaus_kernel);

            if sum(ops.signal_extra_smooth_plot_examples)
                if (n_dset == 1) && (n_pl == 1)
                    figure; 
                    plot(gaus_win*frame_period_sec, gaus_kernel);
                    title(['Gaussian kernel; sigma=' num2str(ops.signal_extra_smooth_sig) 'ms']);
                    xlabel('Time (ms)');
                end
            end

            firing_rate_cut_smooth = conv2(firing_rate_cut, gaus_kernel, 'same');
        else
            firing_rate_cut_smooth = firing_rate_cut;
        end

        if sum(ops.signal_extra_smooth_plot_examples)
            samp_cells = randsample(1:size(firing_rate_cut,1),ops.signal_extra_smooth_plot_examples);
            for n_cell_ind = 1:ops.signal_extra_smooth_plot_examples
                n_cell = samp_cells(n_cell_ind);
                figure; hold on;
                plot(firing_rate_cut(n_cell,:))
                plot(firing_rate_cut_smooth(n_cell,:))
                title(sprintf('Dset %d, plane %d, cell %d', n_dset, n_pl, n_cell));
            end
        end

        % fill in the cut regions
        data.traces_raw{n_dset,n_pl} = if_fill_cuts(traces_raw_cut, cuts_trace);
        data.firing_rate{n_dset,n_pl} = if_fill_cuts(firing_rate_cut, cuts_trace);
        data.firing_rate_smooth{n_dset,n_pl} = if_fill_cuts(firing_rate_cut_smooth, cuts_trace);
        data.num_cells_pl{n_dset,n_pl} = size(data.firing_rate{n_dset,n_pl},1);
        data.stim_frame_index{n_dset,n_pl} = temp_proc_data.stim_frame_index{n_pl};
        cell_plane_indx_pl{n_dset,n_pl} = ones(data.num_cells_pl{n_dset,n_pl},1)*n_pl;
    end
    data.cell_plane_indx{n_dset} = cat(1, cell_plane_indx_pl{:});   
    data.num_cells(n_dset) = sum([data.num_cells_pl{n_dset,:}]);
    if ops.waitbar
        f_waitbar_update(wb, n_dset/num_dsets, sprintf('Preprocessing %d/%d',n_dset,num_dsets));
    end
end
if ops.waitbar
    f_waitbar_close(wb);
end
% check if datasets are equivalent
stim_duration = zeros(num_dsets, 1);
isi = zeros(num_dsets,1);
volume_period = zeros(num_dsets,1);
num_freqs = zeros(num_dsets,1);
for n_dset = 1:num_dsets
    stim_duration(n_dset) = data.proc_data{n_dset,1}.stim_params.stim_duration;
    isi(n_dset) = data.proc_data{n_dset,1}.stim_params.isi;
    volume_period(n_dset) = round(10*data.proc_data{n_dset,1}.frame_data.volume_period_ave)/10;
    num_freqs(n_dset) =  data.proc_data{n_dset,1}.stim_params.num_freqs;
end


if_check_parameter_stability(stim_duration, 'Stim duration');
if_check_parameter_stability(isi, 'isi');
if_check_parameter_stability(volume_period, 'volume_period');
if_check_parameter_stability(num_freqs, 'num_freqs');

end

%%
function firing_rate = if_get_deconvolved_data(temp_data, accepted_cells, ops, frame_period)

traces_raw = temp_data.est.C(accepted_cells,:)+temp_data.est.YrA(accepted_cells,:);
num_cells = size(traces_raw,1);
% so smooth dfdt with default params use for aligning the deconvolution
smooth_dfdt = f_smooth_dfdt3(traces_raw);

if strcmp(ops.signal_inference, 'smooth_dfdt')
    load_params = temp_data.proc.deconv.smooth_dfdt.params;
    params = struct;
    params.do_smooth = load_params.convolve_gaus;
    params.sigma_frames = load_params.gauss_kernel_simga/frame_period;
    params.kernel_window_frames = load_params.gauss_kernel_size;
    params.normalize = load_params.normalize;
    params.rectify = load_params.rectify;
    if load_params.apply_thresh
        params.threshold = load_params.threshold;
    else
        params.threshold = 0;
    end
    firing_rate = f_smooth_dfdt3(traces_raw, params);
elseif strcmp(ops.signal_inference, 'df_f')
    A = temp_data.est.A(:,accepted_cells);
    b = temp_data.est.b;
    C = temp_data.est.C(accepted_cells,:);
    f = temp_data.est.f;
    YrA = temp_data.est.YrA(accepted_cells,:);
    [firing_rate,~] = detrend_df_f_auto_YS(A,b',C,f',YrA);

elseif strcmp(ops.signal_inference, 'c_foopsi') || strcmp(ops.signal_inference, 'MCMC')
    if strcmp(ops.signal_inference, 'c_foopsi')
        firing_rate = temp_data.proc.deconv.c_foopsi.S(accepted_cells,:);
    elseif strcmp(ops.signal_inference, 'MCMC')
        firing_rate = temp_data.proc.deconv.MCMC.S(accepted_cells,:);
    end
    
    % adjust for deconvolution lag
    deconv_lags = zeros(num_cells,1);

    for n_cell = 1:num_cells
        [r,lags] = xcorr(firing_rate{n_cell},smooth_dfdt(n_cell,:),5);
        deconv_lags(n_cell) = lags(r == max(r));
    end

    mean_deconv_lag = round(mean(deconv_lags));
    for n_cell = 1:num_cells
        firing_rate{n_cell} = circshift(firing_rate{n_cell},-1*mean_deconv_lag);
    end
    
    no_deconv_data = false(num_cells,1);
    for n_cell = 1:num_cells
        if isempty(firing_rate{n_cell})
            no_deconv_data(n_cell) = 1;
        end
    end
    if sum(no_deconv_data)
        disp(['Warning: ' cond_name ', dset ' num2str(n_dset) ', plane ' num2str(n_pl) ', ' num2str(sum(no_deconv_data)) '/' num2str(num_cells(n_pl)) ' cells do not have ' ops.signal_inference ' deconvolution data, and will be removed from analysis']);
    end
    firing_rate = cat(1,firing_rate{:});
    firing_rate(no_deconv_data,:) = [];

    
elseif strcmp(ops.signal_inference, 'raw')
    firing_rate = traces_raw;
end

end

function trace_out = if_fill_cuts(trace_in, cuts)
trace_out = zeros(size(trace_in,1), numel(cuts));
trace_out(:,logical(cuts)) = trace_in;
end

function if_check_parameter_stability(parameter, tag)

if ~prod(mode(parameter) == parameter)
    figure;
    plot(parameter, 'o-');
    title(['Warning: ' tag ' varies across datasets'], 'interpreter', 'none');
    ylabel('Stim duration');
    xlabel('Dataset number');
end

end