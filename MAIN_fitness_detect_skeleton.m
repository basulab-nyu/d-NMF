clear;

% Load in Video Y (512 x 512 x t) - Must be double


% Load in Time traces C (k x t) and rois A (512 x 512 x k)


[START_STOP,PEAK] = fitness_detect(C,A,Y);