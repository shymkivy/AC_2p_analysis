function f_dv_update_current_dset_data(app)

n_pl = app.mplSpinner.Value;

num_cells = app.ddata.num_cells_pl{n_pl};
cuts_trace = logical(app.ddata.proc_data{1}.file_cuts_params{n_pl}.vid_cuts_trace);
num_t = numel(cuts_trace);
fr = 1000/double(app.ddata.proc_data{1}.frame_data.volume_period);

%%
accepted_cells = app.ddata.OA_data{n_pl}.proc.comp_accepted;

C = app.ddata.OA_data{n_pl}.est.C;
Yra = app.ddata.OA_data{n_pl}.est.YrA;


if strcmpi(app.DeconvolutionmethodDropDown.Value, 'OA_deconv')
    S = app.ddata.OA_data{n_pl}.est.S;
elseif strcmpi(app.DeconvolutionmethodDropDown.Value, 'MCMC')
    C = app.ddata.OA_data{n_pl}.proc.deconv.MCMC.C;
    S = app.ddata.OA_data{n_pl}.proc.deconv.MCMC.S;
elseif strcmpi(app.DeconvolutionmethodDropDown.Value, 'smooth_dfdt')
    S = app.ddata.OA_data{n_pl}.proc.deconv.smooth_dfdt.S;
end

raw = Yra(accepted_cells,:) + C(accepted_cells,:);

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

%% fill back pulse cuts
raw_full = zeros(num_cells,num_t);
raw_full(:,cuts_trace) = raw;

C_full = zeros(num_cells,num_t);
C_full(:,cuts_trace) = C2;

S_full = zeros(num_cells,num_t);
S_full(:,cuts_trace) = S2;

if app.SmoothCheckBox.Value
    for n_cell = 1:num_cells
        S_full(n_cell,:) = f_smooth_gauss2(S_full(n_cell,:), app.SmoothsigmamsEditField.Value/1000*fr, 0);
    end
end

if app.NormalizespikesCheckBox.Value
    for n_cell = 1:num_cells
        S_full(n_cell,:) = S_full(n_cell,:)/max(S_full(n_cell,:));
    end
end

app.cdata.raw = raw_full;
app.cdata.C = C_full;
app.cdata.S = S_full;

end