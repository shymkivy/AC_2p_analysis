clear;
close all;

f_dir = 'F:\AC_data\caiman_data_missmatch\preprocessing';

flist = dir([f_dir, '\*h5cutsdata.mat']);
f_names = {flist.name};

load_data = load([f_dir, '\', f_names{1}]);


dsall = sum(cat(3,load_data.cuts_data{1}.dsall{:}),3);

figure; plot(dsall)
title(load_data.params.fname, 'interpreter', 'none')









