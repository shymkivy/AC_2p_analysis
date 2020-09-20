function raster_norm = f_normalize(raster, method)

if strcmpi(method, 'norm_full')
    raster_norm = raster - mean(raster,2);
    raster_norm = raster_norm./std(raster_norm,[],2); 
    %firing_rate_cont(isnan(firing_rate_cont)) = 0;
elseif strcmpi(method, 'norm_mean')
    raster_norm = raster - mean(raster,2);
elseif strcmpi(method, 'none')
    raster_norm = raster;
end

end