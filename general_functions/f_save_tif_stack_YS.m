function f_save_tif_stack_YS(stack, name)

% create video average over trials

if nargin < 2 || ~exist('name', 'var')
    name = 'test.tif';
end

% check extension
[~,~,ext] = fileparts(name);
if ~(strcmp(ext, '.tif') || strcmp(ext, '.tiff'))
    name = [name '.tif'];
end

stack = double(stack);

% normalize, range is assumed [0, 1]
stack = stack - min(stack(:));
stack = stack./max(stack(:));

siz = size(stack,3);

hh = waitbar(0, 'Saving Tif stack');
imwrite((stack(:,:,1)), name)
for jj = 2:siz
    imwrite((stack(:,:,jj)), name, 'WriteMode', 'append');
    waitbar(jj/siz, hh);
end
close(hh) 
%implay(norm_diff_vid_stack);

%implay(uint8(norm_diff_vid_stack*(2^8)), 20);

end
