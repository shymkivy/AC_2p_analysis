%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   This scripts takes experiment output files and does preprocessing, by
%   aligning and binning voltage to frames, 
%
%   Input:  XML and CSV files from concatenated videos. the rest files need
%           to contain 'rest' in the name
%           Clictrace outputs
%
%   Output: CSV ile with frame and voltage data (for Jordan's script)
%           _procdata (_info) file with frame data and voltage data
%
%   Preprocessing steps:
%       Extract frame times from XML files
%       Load and process voltage recordings from DAQ
%       Align voltage traces if needed to frames
%       Bin down voltage to frames
%       Concatenate multiple video inputs
%       Locomotion processing
%       Save processed data
%
%   Required functions
%       15x f_s1 function scripts
%       extract_frame_data_from_XML
%       NewCallback_YS
%       align_volt_by_scale_shift2
%       ask_yes_no_fig
%       
%   Last update: Yuriy 12/2/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script
clear;
close all;

%% file path
addpath([ pwd '\s1_functions'])
addpath([ pwd '\general_functions'])

ops.exp_date = '10_14_19';
%ops.file_core=['DAM_ammn_2_dplanes1_'  ops.exp_date ''];
ops.file_core=['DF_ammn2_5plt_2plx_'  ops.exp_date ''];
%ops.file_dir = 'L:\data\Auditory\caiman_out_multiplane';
ops.file_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data\preprocessing';
%ops.exp_date = '';
%ops.file_core=['rfp-asap3-grat-008'  ops.exp_date ''];
%ops.file_dir = 'C:\Users\ys2605\Desktop\victor crap data\with_pulses\rfp-asap3-grat-008';


% voltage csv and xml and traces raw base
ops.files_volt_in = {['raw_' ops.file_core]};
% files to be concetentated
%files_in = {['audio_grating1raw_' exp_date]};
              %['rest7_' exp_date];
              %['rest6_' exp_date];
              %['grating1_' exp_date];
              %['vmmn2_' exp_date]};

%% params if using on acid vedo cuts file
ops.ca_processing = 'onacid';          % options: 'onacid', 'clicktrace', 'raw_movie';

% for multiplane dataset set number of planes, otherwise 1;
ops.num_planes = 5;   

% plot extra details along the way
ops.plot_details = 1;

% What to preprocess
ops.processing_type = 2;
% = 1; rest vs grating
% = 2; control, MMN, flipMMN (uses stim chan for getting stim_times)
% = 3; auditory freq grating; (needs TDT_volt for stim_times)

% DAQ voltage channels recording order
ops.parameters.stimchan = 1;
ops.parameters.ledchan = 2;
ops.parameters.movchan = 3;
ops.parameters.TDT_volt_chan = 4;

% which voltage channel to use for alignment? 1 for video 2 for auditory
ops.align_to_channel = 2;

ops.alignment_method = 'peak_onsets_scale_only'; 
% options: 'xcorr', 'peak_onsets', 'peak_onsets_shift_only', 'peak_onsets_scale_only', 'manual'
% options: 


ops.exp_window_selection = 'auto';  % options: 'auto', 'manual'

% voltage output for jordan
ops.bin_csv_out = 0;

ops.auto_loco_thresh = 0; % 0 for manual, >0 for auto;

%% process ops
ops = f_s1_process_ops(ops);

%% Load presaved parameters if exist
data = struct;
if exist(ops.file_save_path_full_processing_params, 'file')
    temp_data = load(ops.file_save_path_full_processing_params);
    data_fieldnames = fieldnames(temp_data);
    for n_fl = 1:numel(data_fieldnames)
        data.(data_fieldnames{n_fl}) = temp_data.(data_fieldnames{n_fl});
    end
    clear n_fl temp_data data_fieldnames;
end

%% Load data from from XML file
data = f_s1_load_frame_times(data, ops);

%% Load and process voltage traces
data = f_s1_load_voltage_data(data, ops);

%% Load calcium data
data = f_s1_load_ca_data(data, ops);

%% Load stim data (matlab output)
[data, ops] = f_s1_load_stim_params(data, ops);

%% Get alignment information
data = f_s1_get_alignment_info(data, ops);

%% Do the alignment, trimming, and concatenation
data = f_s1_apply_alignment(data, ops);

%% for auditory script generate stim start times in ms
data = f_s1_process_volt(data, ops);

%% Define experimental windows
data = f_s1_define_exp_windows(data, ops);

%% Stim times calculatoin
data = f_s1_process_stim_volt(data, ops);

%% quality plots
% plot some things to check alignment and preprocessing
f_s1_plot_processing_quality(data, ops);

%% save binned down voltage data
fprintf('Saving files:\n');

save([ops.file_save_path_full '_processed_data.mat'], 'data', 'ops', '-v7.3');

fprintf('Saved %s_processed_data.mat\n', ops.file_save_path_full);

if ops.bin_csv_out
    for n_pl = 1:ops.num_planes
        dataall = volt_data{n_pl}';
        save([ops.file_path_full_ca{n_pl} '_bin.csv'],'dataall','-ascii');
        fprintf('%s.csv\n', ops.file_save_path_full);

        if ops.OnACID
            dataall_cut = dataall(logical(data.file_cuts_params{n_pl}.vid_cuts_trace));
            save([ops.file_path_full_ca{n_pl} '_bin_cut.csv'],'dataall_cut','-ascii');
            clear dataall_cut;
        end
    end
    clear dataall;
end

%  print errors if present
if ~isempty(ops.errors)
    disp('Erorrs:');
    for ii = 1:numel(ops.errors)
        disp(ops.errors{ii});
    end
end

fprintf('Done\n');

