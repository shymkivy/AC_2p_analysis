function Y = f_collect_prairie_tiffs2(load_stack, key)
% load the PrairieView output tif images and combine them into single
% movie stack h5 or tif. h5 is much faster

num_cat = numel(load_stack);

T = zeros(numel(load_stack),1);
dir_names = cell(numel(load_stack),1);

% check if all prairie directories exists

for n_cat = 1:num_cat
    if ~exist(load_stack{n_cat}, 'dir')
        Error(['Pipeline warning: Directory ' load_stack{n_cat} ' does not exist']);
        continue;
    end

    dir_names{n_cat} = dir([load_stack{n_cat}, '\*.tif']);
    
    if exist('key', 'var')
        read_ind = false(numel(dir_names{n_cat}),1);
        for ii = 1:numel(dir_names{n_cat})
            read_ind(ii) = ~isempty(strfind(dir_names{n_cat}(ii).name,key));
        end
        T(n_cat) = sum(read_ind);
        dir_names{n_cat}(read_ind) = [];
    else
        T(n_cat) = numel(dir_names{n_cat});
    end
    
    
end

tremp_frame = imread([load_stack{1}, '\',  dir_names{1}(1).name], 'tif');
dim = size(tremp_frame);
T_tot = sum(T);
Y = zeros([dim, T_tot], 'uint16');

if T_tot == 0
    Error('Pipeline Error: No tiff files in specified folters');
end

shift = 0;
for n_cat = 1:num_cat

    hh = waitbar(0, ['Loading Prairie tiffs; file ', num2str(n_cat), ' of ', num2str(num_cat)]);
    for n_frame = 1:T(n_cat)
        Y(:,:,n_frame+shift) = imread([load_stack{n_cat}, '\', dir_names{n_cat}(n_frame).name], 'tif');
        waitbar((n_frame+shift)/T_tot, hh);
    end
    shift = shift+T(n_cat);
    close(hh)
end

end