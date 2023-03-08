function [options] = defaultOptions_mcbDNMF()
    %[options] = defaultOptions_mcbDNMF()
    % Size of image patches
    options.patchSize = [64 64];        % Size of image patches
    options.stride = 56;                % Offset of image patches

    % Initial ROI detection parameters
    options.DETREND_FRAMES = 45;        % Number of frames to compute baseline F0 over
    options.med_opt = false;             % If true, baseline F0 is computed as median; otherwise computed as minimum
    options.thr = 3;                    % Threshold for active pixels 
    options.filtSize = 3;               % Median filter size for initial ROI detection
    options.overlapThr = 0.5;           % Spatial overlap merge threshold 
    options.sizeRange = [30 2000];      % Allowable size range of valid ROIs
    options.eta = 0.01;                 % Temporal regularization weight
    options.beta = 0.5;                 % Spatial regularization weight
    options.med_opt = false;            % Option to compute baseline fluoresence as median rather than minimum. Recommended for non-downsampled data

    % ROI cleanup parameters
    options.thr_method = 'quant';       % Method of thresholding ROIs: 'max' or 'quant'
    options.quantileThr = 0.9;          % Quantile threshold for thresholding using 'quant'
    options.maxthr = 0.2;               % Max threshold for thresholding using 'max'
    options.final_C = true;             % Whether or not to recompute traces C after ROI cleanup

    % ROI validity & merging parameters
    options.minSkew = 0;                % Minimum skew of temporal trace of valid ROIs
    options.shapeThr = 0.5;             % Correlation threshold for ROIs
    options.temporalCorrThr = 0.9;     % Temporal correlation merge threshold 

end