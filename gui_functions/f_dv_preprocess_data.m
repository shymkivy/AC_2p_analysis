function data = f_dv_preprocess_data(data, ops)

disp('Preprocessing...');

num_dsets = numel(data.area);
if ops.waitbar
    if isfield(ops, 'app')
        app = ops.app;
    else
        app = [];
    end
    wb = f_waitbar_initialize(app, 'Preprocessing data...');
end
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    proc_data = ddata.proc_data{1};
    %frame_period_sec = proc_data.frame_data.volume_period/1000;
    shift_stim = 0;
    if isfield(proc_data, 'trial_types')
        trial_types = proc_data.trial_types;
        if sum(sum(unique(trial_types) == 1:10,2)) == 8
            idx1 = trial_types == 8;
            trial_types(idx1) = trial_types(idx1)+2;
            idx1 = logical(sum(trial_types == 4:7,2));
            trial_types(idx1) = trial_types(idx1)+1;
            shift_stim = 1;
        end
        data.trial_types{n_dset} = trial_types;
    end

    if isfield(proc_data, 'stim_params')
        if isfield(proc_data.stim_params, 'MMN_freq')
            mmn_freq = proc_data.stim_params.MMN_freq;
        elseif isfield(proc_data, 'MMN_orientations')
            mmn_freq = proc_data.MMN_orientations;
        else
            mmn_freq = [];
            sprintf('No mmn freq in %s', ddata.dset_name_full{1});
        end
        if shift_stim
            idx1 = mmn_freq == 8;
            mmn_freq(idx1) = mmn_freq(idx1)+2;
            idx1 = logical(sum(mmn_freq' == 4:7,2));
            mmn_freq(idx1) = mmn_freq(idx1)+1;
        end
        data.MMN_freq{n_dset} = mmn_freq;
        
    end
    
    cell_plane_indx_pl = cell(data.num_planes(n_dset),1);
    disp(ddata.dset_name_full{1})
    % pull out data
    for n_pl = 1:data.num_planes(n_dset)
        temp_OA_data = data.OA_data{n_dset,n_pl};
        proc_data = data.proc_data{n_dset};
        % extra cSNR threshold
        SNR_accept = temp_OA_data.proc.SNR2_vals >= ops.extra_SNR_thresh;
        accept_cell = and(SNR_accept,temp_OA_data.proc.comp_accepted);           
        traces_raw_cut = temp_OA_data.est.C(accept_cell,:)+temp_OA_data.est.YrA(accept_cell,:);

        % fill in the cut regions
        cuts_trace = proc_data.file_cuts_params{n_pl}.vid_cuts_trace;
        data.traces_raw{n_dset,n_pl} = if_fill_cuts(traces_raw_cut, cuts_trace);
        %data.firing_rate{n_dset,n_pl} = if_fill_cuts(firing_rate_cut, cuts_trace);
        %data.firing_rate_smooth{n_dset,n_pl} = if_fill_cuts(firing_rate_cut_smooth, cuts_trace);
        data.num_cells_pl{n_dset,n_pl} = size(data.traces_raw{n_dset,n_pl},1);
        if isfield(proc_data, 'stim_times_frame')
            data.stim_frame_index{n_dset,n_pl} = proc_data.stim_times_frame{1, n_pl};
        elseif isfield(proc_data, 'stim_frame_index')
            data.stim_frame_index{n_dset,n_pl} = proc_data.stim_frame_index{n_pl};
        end
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
volume_period = zeros(num_dsets,1);
for n_dset = 1:num_dsets
    volume_period(n_dset) = round(10*data.proc_data{n_dset,1}.frame_data.volume_period_ave)/10;
end
if_check_parameter_stability(volume_period, 'volume_period');

stim_duration = zeros(num_dsets, 1);
isi = zeros(num_dsets,1);
num_freqs = zeros(num_dsets,1);
has_data = false(num_dsets,1);
for n_dset = 1:num_dsets
    if isfield(data.proc_data{n_dset,1}, 'stim_params')
        has_data(n_dset) = 1;
        stim_duration(n_dset) = data.proc_data{n_dset,1}.stim_params.stim_duration;
        isi(n_dset) = data.proc_data{n_dset,1}.stim_params.isi;
        num_freqs(n_dset) =  data.proc_data{n_dset,1}.stim_params.num_freqs;
    end
end

if_check_parameter_stability(stim_duration(has_data), 'Stim duration');
if_check_parameter_stability(isi(has_data), 'isi');
if_check_parameter_stability(num_freqs(has_data), 'num_freqs');


end

%%
function firing_rate = if_get_deconvolved_data(temp_data, accepted_cells, ops, frame_period)

traces_raw = temp_data.est.C(accepted_cells,:)+temp_data.est.YrA(accepted_cells,:);
num_cells = size(traces_raw,1);
% so smooth dfdt with default params use for aligning the deconvolution
smooth_dfdt = f_smooth_dfdt3(traces_raw);

if strcmpi(ops.signal_inference, 'smooth_dfdt')
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
elseif strcmpi(ops.signal_inference, 'df_f')
    A = temp_data.est.A(:,accepted_cells);
    b = temp_data.est.b;
    C = temp_data.est.C(accepted_cells,:);
    f = temp_data.est.f;
    YrA = temp_data.est.YrA(accepted_cells,:);
    [firing_rate,~] = detrend_df_f_auto_YS(A,b',C,f',YrA);

elseif strcmpi(ops.signal_inference, 'c_foopsi') || strcmpi(ops.signal_inference, 'MCMC')
    if strcmpi(ops.signal_inference, 'c_foopsi')
        firing_rate = temp_data.proc.deconv.c_foopsi.S(accepted_cells,:);
    elseif strcmpi(ops.signal_inference, 'MCMC')
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

    
elseif strcmpi(ops.signal_inference, 'raw')
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
    ylabel(tag, 'interpreter', 'none');
    xlabel('Dataset number');
end

end