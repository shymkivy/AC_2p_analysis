function imData=bigread(path_to_file,varargin)
%reads tiff files in Matlab bigger than 4GB.  
%hacked and extended version of a partial solution posted on Matlab Central (http://www.mathworks.com/matlabcentral/answers/108021-matlab-only-opens-first-frame-of-multi-page-tiff-stack)
%Darcy Peterka 2014, v1.0
%Program checks for bit depth, whether int or float, and byte order.  Assumes uncompressed, non-negative (i.e. unsigned) data.
%
% Usage:  my_data=bigread('path_to_data_file;);
% "my_data" will be your [M,N,frames] array.
if ~isempty(varargin)
    ml=varargin(1);
    maxframes=ml{1};
    disp(maxframes)
else
    maxframes=9999999;
end
%get image info
info = imfinfo(path_to_file);
he=info.ImageDescription;
numFramesStr = regexp(he, 'images=(\d*)', 'tokens');
if isempty(numFramesStr)
    blah=size(info);
    numFrames= blah(1);
else
    numFrames = str2double(numFramesStr{1}{1});
end

bd=info.BitDepth;
he=info.ByteOrder;
bo=strcmp(he,'big-endian');
if (bd==32)
	form='single';
elseif (bd==16)
		form='uint16';
elseif (bd==8)
		form='uint8';
end

% Use low-level File I/O to read the file
fp = fopen(path_to_file , 'rb');
% The StripOffsets field provides the offset to the first strip. Based on
% the INFO for this file, each image consists of 1 strip.
he=info(1).StripOffsets;
% if (max(size(he))>1)
%     error('Tiff has non-uniform offsets - I''m not coded to handle this...sorry');
% end
fseek(fp, he, 'bof');
framenum=numFrames;
imData=cell(1,framenum);
he_w=info.Width;
he_h=info.Height;
if framenum>maxframes
    framenum1=maxframes;
else
    framenum1=framenum;
end
h = waitbar(0, 'Loading Video...');
if(bo)
	for cnt = 1:framenum1
        %cnt
		tmp1 = fread(fp, [he_w he_h], form, 0, 'ieee-be')';
        imData{cnt}=cast(tmp1,'uint16');
        waitbar(cnt/framenum1);
	end
else
	for cnt = 1:framenum1
         %cnt
		tmp1 = fread(fp, [he_w he_h], form, 0, 'ieee-le')';
        imData{cnt}=cast(tmp1,'uint16');
        waitbar(cnt/framenum1);
	end
end
close(h);
imData=cell2mat(imData);
imData=reshape(imData,[he_h,he_w,framenum1]);
fclose(fp);