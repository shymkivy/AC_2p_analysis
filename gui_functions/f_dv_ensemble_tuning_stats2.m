function ens_tuning_stats = f_dv_ensemble_tuning_stats2(app, params)

ens_out = params.ensembles.ens_out;
num_ens = ens_out.num_comps;

if isfield(params, 'ensemble_stats')
    accepted_ens = params.ensemble_stats.accepted_ensembles;
else
    accepted_ens = ones(num_ens, 1);
end

ens_scores = ens_out.scores(accepted_ens,:);

cdata.S_sm = ens_scores;
cdata.volume_period = params.cdata(1,:).volume_period;
cdata.num_cells = size(ens_scores,1);
cdata.accepted_cells = accepted_ens;

ens_params = params;
ens_params.cdata = cdata;

ens_tuning_stats = f_dv_compute_stats_core(app, ens_params);

end