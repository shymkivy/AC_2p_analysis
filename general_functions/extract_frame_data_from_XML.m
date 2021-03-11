function [data_out, errors] = extract_frame_data_from_XML(file_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               
%   Purpose: This function extracts frame data from XML file    
%   Date: 3/8/18                                                
%                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


fprintf('Loading XML file %s...\n', file_name);

% collect errors in this
errors = {};

[path1,fname1,~] = fileparts(file_name);
file_struct = xml2struct([path1, '\', fname1, '.xml']);

disp('Processing XML file...');

% extract useful data
data_out.xml_version = file_struct.PVScan.Attributes.version;
data_out.frame_period = 1000*str2double(file_struct.PVScan.Sequence.Frame{2}.Attributes.relativeTime);
data_out.volt_start = 1000*str2double(file_struct.PVScan.Sequence.VoltageRecording.Attributes.absoluteTime);
data_out.frame_start = 1000*str2double(file_struct.PVScan.Sequence.Frame{1}.Attributes.absoluteTime);

% counts number of frames
data_out.total_frames = length(file_struct.PVScan.Sequence.Frame);

% extract frame times
data_out.frame_times = zeros(data_out.total_frames, 1);
for jj = 1:data_out.total_frames
    data_out.frame_times(jj) = 1000*str2double(file_struct.PVScan.Sequence.Frame{jj}.Attributes.relativeTime);
end

% Fix some specific bugs belowbugs below

% fix bug if version 5.4 alpha, prame times are double with 256 2ave scan
try
    if sum(data_out.xml_version == '5.4.64.63') == 9
        if data_out.frame_period > 50 % check if larger than 50 ms
            data_out.frame_period = data_out.frame_period/2;
            data_out.frame_times = data_out.frame_times/2;
            errors = [errors; {sprintf('Fixed 5.4alpha bug for %s',file_name)}];
        end
    end
end

clear file_struct;
end