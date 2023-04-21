clear;
tiff_files = {   'C:\Users\Jason\Documents\MATLAB\DNMF_Review-main\Data\61\61D5._Tsub_mean.tif'};
dnmf_files = {   'C:\Users\Jason\Documents\MATLAB\DNMF_Review-main\Data\61\DNMF_Out.mat'};

% Load in Video Y (512 x 512 x t) - Must be double
Y = double(bigread2(tiff_files{1}));

% Load in Time traces C (k x t) and rois A (512 x 512 x k)
load(dnmf_files{1});
C = Cs;
A = reshape(full(cROIs),[512 512 size(cROIs,2)]);
[START_STOP,PEAK] = fitness_detect(C,A,Y);