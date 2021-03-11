function data_out = extract_frame_data_from_XML2(file_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               
%   Purpose: This function extracts frame data from XML file, faster than
%   previous version
%
%   Date: 11/25/19            
%   Yuriy
%                                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf('Loading XML file %s...\n', file_name);
    
    % collect errors in this
errors = {};


[path1,fname1,~] = fileparts(file_name);

dataXML = xmlread([path1, '\', fname1, '.xml']);


data_out.xml_version = string(dataXML.getElementsByTagName('PVScan').item(0).getAttribute('version'));
data_out.volt_start = 1000*str2double(dataXML.getElementsByTagName('VoltageRecording').item(0).getAttribute('absoluteTime'));
data_out.frame_start = 1000*str2double(dataXML.getElementsByTagName('Frame').item(0).getAttribute('absoluteTime'));
data_out.trigger_mode = string(dataXML.getElementsByTagName('VoltageRecording').item(0).getAttribute('triggerMode'));

AllFramesList = dataXML.getElementsByTagName('Frame');

data_out.total_frames = AllFramesList.getLength;
data_out.frame_times = zeros(data_out.total_frames, 1);
for n_fr = 1:data_out.total_frames
    data_out.frame_times(n_fr) = 1000*str2double(AllFramesList.item(n_fr-1).getAttribute('relativeTime'));
end

data_out.frame_period = mean(diff(data_out.frame_times));

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

data_out.errors = errors;

end