function cdata = f_dv_compute_cdata(app, params)

n_pl = params.n_pl;
n_dset = params.n_dset;

ddata = app.data(n_dset,:);

cuts_trace = logical(ddata.proc_data{1}.file_cuts_params{n_pl}.vid_cuts_trace);
num_t = numel(cuts_trace);
fr = 1000/double(ddata.proc_data{1}.frame_data.volume_period);

%%
accepted_cells = ddata.OA_data{n_pl}.proc.comp_accepted;

% filter out cells with missing sig
sig_frac = 1 - ddata.OA_data{n_pl}.proc.num_zeros/ddata.OA_data{n_pl}.proc.num_frames;
accepted_cells(sig_frac<.9) = 0;

num_cells = sum(accepted_cells);

%%
C = ddata.OA_data{n_pl}.est.C;
Yra = ddata.OA_data{n_pl}.est.YrA;
raw = Yra + C;
raw2 = raw(accepted_cells,:);

if strcmpi(params.deconvolution, 'OA_deconv')
    S = ddata.OA_data{n_pl}.est.S;
elseif strcmpi(params.deconvolution, 'MCMC')
    C = ddata.OA_data{n_pl}.proc.deconv.MCMC.C;
    S = ddata.OA_data{n_pl}.proc.deconv.MCMC.S;
elseif strcmpi(params.deconvolution, 'smooth_dfdt')
    sigma1 = ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.params.gauss_kernel_simga;
    sigma_frames = sigma1/1000*fr;
    do_smooth = 1;
    normalize1 = 0;
    rectify1 = params.rectify_spikes;
    S = f_smooth_dfdt3(raw, do_smooth, sigma_frames, normalize1, rectify1);
    %S = app.ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.S;
end

if iscell(C(1,:))
    C_temp = C(accepted_cells);
    C2 = cat(1,C_temp{:});
else
    C2 = C(accepted_cells,:);
end

if iscell(S(1,:))
    S_temp = S(accepted_cells);
    S2 = cat(1,S_temp{:});
else
    S2 = S(accepted_cells,:);
end

if params.subtract_mean_spikes
    S2 = S2 - mean(S2,2);
end

if params.normalize_max_spikes
    S2 = S2./max(S2,[],2);
end

%% fill back pulse cuts
raw_full = zeros(num_cells,num_t);
raw_full(:,cuts_trace) = raw2;

C_full = zeros(num_cells,num_t);
C_full(:,cuts_trace) = C2;

S_full = zeros(num_cells,num_t);
S_full(:,cuts_trace) = S2;

if params.smooth
    for n_cell = 1:num_cells
        S_full(n_cell,:) = f_smooth_gauss2(S_full(n_cell,:), params.smooth_sigma/1000*fr, 0);
    end
end

cdata.raw = raw_full;
cdata.C = C_full;
cdata.S = S_full;
cdata.accepted_cells = accepted_cells;
cdata.num_cells = num_cells;
cdata.volume_period = ddata.proc_data{1}.frame_data.volume_period;

end