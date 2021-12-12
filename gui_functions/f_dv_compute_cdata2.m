function cdata = f_dv_compute_cdata2(app, params)

n_dset = params.n_dset;
ddata = app.data(n_dset,:);

data_gr = app.SelectdatagroupButtonGroup.SelectedObject.Text;
if strcmpi(data_gr, 'plane')
    planes1 = params.n_pl;
else
    planes1 = 1:ddata.num_planes;
end

fr = 1000/double(ddata.proc_data{1}.frame_data.volume_period);

%%
num_planes = numel(planes1);

accepted_cells = cell(num_planes,1);
num_cells = zeros(num_planes,1);
C = cell(num_planes,1);
Yra = cell(num_planes,1);
raw = cell(num_planes,1);
S = cell(num_planes,1);


for n_pl = planes1
    accepted_cells2 = ddata.OA_data{n_pl}.proc.comp_accepted;
    sig_frac = 1 - ddata.OA_data{n_pl}.proc.num_zeros/ddata.OA_data{n_pl}.proc.num_frames;
    accepted_cells2(sig_frac<.9) = 0;
    num_cells2 = sum(accepted_cells2);
    
    C2 = ddata.OA_data{n_pl}.est.C;
    Yra2 = ddata.OA_data{n_pl}.est.YrA;
    raw2 = Yra2 + C2;
    raw3 = raw2(accepted_cells2,:);

    if strcmpi(params.deconvolution, 'OA_deconv')
        S2 = ddata.OA_data{n_pl}.est.S;
    elseif strcmpi(params.deconvolution, 'MCMC')
        C2 = ddata.OA_data{n_pl}.proc.deconv.MCMC.C;
        S2 = ddata.OA_data{n_pl}.proc.deconv.MCMC.S;
    elseif strcmpi(params.deconvolution, 'smooth_dfdt')
        sigma1 = ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.params.gauss_kernel_simga;
        sigma_frames = sigma1/1000*fr;
        do_smooth = 1;
        normalize1 = 0;
        rectify1 = params.rectify_spikes;
        S2 = f_smooth_dfdt3(raw2, do_smooth, sigma_frames, normalize1, rectify1);
        %S = app.ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.S;
    end
    
    if iscell(C2(1,:))
        C_temp = C2(accepted_cells2);
        C3 = cat(1,C_temp{:});
    else
        C3 = C2(accepted_cells2,:);
    end
    
    if iscell(S2(1,:))
        S_temp = S2(accepted_cells2);
        S3 = cat(1,S_temp{:});
    else
        S3 = S2(accepted_cells2,:);
    end
    
    if params.subtract_mean_spikes
        S3 = S3 - mean(S3,2);
    end

    if params.normalize_max_spikes
        S3 = S3./max(S3,[],2);
    end

    
    cuts_trace = logical(ddata.proc_data{1}.file_cuts_params{n_pl}.vid_cuts_trace);
    num_t = numel(cuts_trace);
    
    accepted_cells{n_pl} = accepted_cells2;
    num_cells(n_pl) = num_cells2;
    
    C{n_pl} = zeros(num_cells2,num_t);
    C{n_pl}(:,cuts_trace) = C3;
    
    Yra{n_pl} = zeros(num_cells2,num_t);
    Yra{n_pl}(:,cuts_trace) = Yra2(accepted_cells2,:);
    
    raw{n_pl} = zeros(num_cells2,num_t);
    raw{n_pl}(:,cuts_trace) = raw3;
    
    S{n_pl} = zeros(num_cells2,num_t);
    S{n_pl}(:,cuts_trace) = S3;
    
    if params.smooth
        for n_cell = 1:num_cells2
            S{n_pl}(n_cell,:) = f_smooth_gauss2(S{n_pl}(n_cell,:), params.smooth_sigma/1000*fr, 0);
        end
    end
    
end


%%

cdata.raw = raw;
cdata.C = C;
cdata.S = S;
cdata.accepted_cells = accepted_cells;
cdata.num_cells = num_cells;
cdata.volume_period = ddata.proc_data{1}.frame_data.volume_period;

end