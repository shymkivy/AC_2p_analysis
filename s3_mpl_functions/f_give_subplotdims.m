function [m,n] = f_give_subplotdims(num_plots)

if num_plots < 4
    m = 1;
    n = num_plots;
elseif num_plots < 11
    m = 2;
    n = ceil(num_plots/2);
elseif num_plots < 11
    m = 3;
    n = ceil(num_plots/3);
else
    m = ceil(num_plots/5);
    n = 5;
end


end