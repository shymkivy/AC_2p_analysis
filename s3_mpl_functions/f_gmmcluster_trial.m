function f_gmmcluster_trial(X, params)


method = params.method;
metric = params.metric;
num_cells = size(X,1);

%% estimate num clust
if ~isfield(params, 'num_clust')
    
    if num_cells > 5
        k_list = 1:6;
    else
        k_list = 1:size(X,1);
    end

    E = evalclusters(X,'gmdistribution','silhouette','klist',k_list, 'Distance', metric);
    num_clust = E.OptimalK;
    
else
    num_clust = params.num_clust;
end









end