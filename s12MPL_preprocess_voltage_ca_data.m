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

ops.exp_date = '1_31_22';
%ops.file_core=['DAM_ammn_2_dplanes1_'  ops.exp_date ''];
ops.file_core=['M101_im1_AC_ammn1_'  ops.exp_date ''];
%ops.file_dir = 'L:\data\Auditory\caiman_out_multiplane';
ops.file_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\\caiman_data_dream\preprocessing';
%ops.exp_date = '';
%ops.file_core=['rfp-asap3-grat-008'  ops.exp_date ''];
%ops.file_dir = 'C:\Users\ys2605\Desktop\victor crap data\with_pulses\rfp-asap3-grat-008';


% voltage csv and xml and traces raw base
%ops.files_volt_in = {['raw_' ops.file_core]};
ops.files_volt_in = {['' ops.file_core '_prairie']};
% files to be concetentated
%files_in = {['audio_grating1raw_' exp_date]};
              %['rest7_' exp_date];
              %['rest6_' exp_date];
              %['grating1_' exp_date];
              %['vmmn2_' exp_date]};

%%
f_s12_preprocess_voltage_ca_data(ops);
    
%%
fprintf('Done\n');

