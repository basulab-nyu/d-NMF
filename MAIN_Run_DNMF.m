clear;
addpath(genpath('.'));
%% Set file path
    files{1} = {'C:\Users\Jason\Documents\MATLAB\Dugtrio\2020-08-20-002\MC_Video_TSub_nonrigid.tif'};

    
    
%% Set options - Downsampled (1.5 Hz)
options.maxVal = 2^13;              % Maximum pixel intensity. Set anything above to 0


% Size of image patches
options.patchSize = [64 64];        % Size of image patches
options.stride = 56;                % Offset of image patches
options.loadInPatches = false;      % Option to load in data in patches (true) or all at once (false).Set to true to reduce memory requirements at the expense of speed

% Initial ROI detection parameters
options.spatial_smooth = 0;         % Sigma of spatial gaussian smoothing kernel, in pixels
options.temporal_filter_corner = 0.5;   % Corner cutoff frequency for temporal low-pass filter, in Hz
options.fs = 1.5;                   % Sampling rate (Hz)
options.SNR_thr = 8;                % Threshold on signal-to-noise ratio for screening initial ROIs
options.DETREND_FRAMES = round(30*options.fs);        % Number of frames to compute baseline F0 over
options.thr = 3;                    % Threshold for active pixels 
options.filtSize = 1;               % Median filter size for initial ROI detection
options.overlapThr = 0.4;           % Spatial overlap merge threshold 
options.sizeRange = [30 2000];      % Allowable size range of valid ROIs
options.eta = 10;                   % Temporal regularization weight
options.beta = 5e5;                 % Spatial regularization weight

options.DOWNSAMPLE_FACTOR = 1;      % Amount to temporally downsample movie for initialization
options.med_opt = false;            % Option to compute baseline fluoresence as median rather than minimum. Recommended true for non-downsampled data

% ROI cleanup parameters
options.thr_method = 'max';         % Method of thresholding ROIs: 'max' or 'quant'
options.quantileThr = 0.9;          % Quantile threshold for thresholding using 'quant'
options.maxthr = 0.2;               % Max threshold for thresholding using 'max'
options.final_C = true;             % Whether or not to recompute traces C after ROI cleanup

% ROI validity & merging parameters
options.minSkew = 0;                % Minimum skew of temporal trace of valid ROIs
options.shapeThr = 0.2;             % Correlation threshold for ROIs



%% Set options - Full Rate (30 Hz)
options.maxVal = 2^13;              % Maximum pixel intensity. Set anything above to 0

% Size of image patches
options.patchSize = [64 64];        % Size of image patches
options.stride = 56;                % Offset of image patches
options.loadInPatches = true;       % Option to load in data in patches (true) or all at once (false).Set to true to reduce memory requirements at the expense of speed

% Initial ROI detection parameters
options.spatial_smooth = 0;         % Sigma of spatial gaussian smoothing kernel, in pixels
options.temporal_filter_corner = 0; % Corner cutoff frequency for temporal low-pass filter, in Hz
options.fs = 30;                    % Sampling rate (Hz)
options.SNR_thr = 8;                % Threshold on signal-to-noise ratio for screening initial ROIs
options.DETREND_FRAMES = round(30*options.fs);        % Number of frames to compute baseline F0 over
options.thr = 3;                    % Threshold for active pixels 
options.filtSize = 1;               % Median filter size for initial ROI detection
options.overlapThr = 0.4;           % Spatial overlap merge threshold 
options.sizeRange = [30 2000];      % Allowable size range of valid ROIs
options.eta = 10;                 % Temporal regularization weight
options.beta = 1e6;                 % Spatial regularization weight

options.DOWNSAMPLE_FACTOR = 20;     % Amount to temporally downsample movie for initialization
options.med_opt = false;            % Option to compute baseline fluoresence as median rather than minimum. Recommended true for non-downsampled data
% For this parameter set, med_opt should be false, as we compute baseline
% fluorescence on downsampled data, so we use the minimum. In general, for
% high-rate data (>5 Hz) med_opt is recommended to be true.


% ROI cleanup parameters
options.thr_method = 'max';         % Method of thresholding ROIs: 'max' or 'quant'
options.quantileThr = 0.9;          % Quantile threshold for thresholding using 'quant'
options.maxthr = 0.2;               % Max threshold for thresholding using 'max'
options.final_C = true;             % Whether or not to recompute traces C after ROI cleanup

% ROI validity & merging parameters
options.minSkew = 0;                % Minimum skew of temporal trace of valid ROIs
options.shapeThr = 0.2;             % Correlation threshold for ROIs

%% Main Loop
tempFolder = '.';                   % Temporary folder to write ROIs to before zipping them up to target folder
for i_file = length(files)
    thisFile = files{i_file};
    
    [cROIs, Cs, coherence, skew, sz, patchID, tElapsed, dimensions] = mcb_DNMF(thisFile, options);
    
    folder = fileparts(thisFile);
    outFolder = folder;
    if(~exist(outFolder,'dir'))
        mkdir(outFolder);
    end
    fprintf('Saving %d ROIs...', size(cROIs, 2));
    save(fullfile(outFolder, 'DNMF_Out.mat'), 'cROIs', 'Cs', 'coherence', 'skew', 'sz', 'patchID', 'options', 'tElapsed', '-v7.3');  
        
    if(~isempty(tempFolder))
        outFolder_original = outFolder;
        outFolder = tempFolder;
    end
    outputFolder = fullfile(outFolder, 'ROI_Set\');
    if(~exist(outputFolder, 'dir'))
        mkdir(outputFolder);
    end
    [~,max_idx] = max(Cs, [], 2);
    for i_roi = 1:size(cROIs, 2)   
        temp = reshape(full(cROIs(:, i_roi)), [512 512]);
        [a, b] = find(imdilate(temp>0, ones(3)));
        c = boundary(a, b, 0.95);

        writeImageJROI_3([a(c) b(c)], 4, max_idx(i_roi), sprintf('r%04d',i_roi), outputFolder);
    end

    zip(strrep(outputFolder, 'ROI_Set\', 'RoiSet_DNMF'), '*', outputFolder);
    rmdir(outputFolder, 's');
    if(~isempty(tempFolder))
        movefile(fullfile(outFolder, 'RoiSet_DNMF.zip'), fullfile(outFolder_original, 'RoiSet_Auto.zip'));
    end
end