function f_save_tif_stack2_YS(stack, filename)
% updated 9/2020

% check extension
[filepath,name,ext] = fileparts(filename);
if ~(strcmpi(ext, '.tif') || strcmpi(ext, '.tiff'))
    ext = '.tif';
    filename = [filepath '\' name ext];
end

data_type = class(stack);

if ~strcmpi(data_type(1:min(4,numel(data_type))),'uint')
    stack = double(stack);
    % normalize, range is assumed [0, 1]
    stack = stack - min(stack(:));
    stack = stack./max(stack(:));
end

T = size(stack,3);

num_err = 0;

hh = waitbar(0, 'Saving Tif stack');
imwrite((stack(:,:,1)), filename, 'tif');
for jj = 2:T
    n_try = 1;
    while n_try<10
        try
            imwrite((stack(:,:,jj)), filename, 'tif', 'WriteMode', 'append');
            n_try = 10;
        catch
            n_try = n_try + 1;
            num_err = num_err + 1;
        end
    end
    waitbar(jj/T, hh, sprintf('Saving Tif stack; error count = %d', num_err));
end
close(hh) 
%implay(norm_diff_vid_stack);

%implay(uint8(norm_diff_vid_stack*(2^8)), 20);

end
