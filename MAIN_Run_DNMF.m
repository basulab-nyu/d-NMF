clear;
addpath(genpath('.'));
%% Set file path
files = {   'H:\MATLAB\Moore_2023_Code\DNMF_Alpha\Data\495\D02\495D2G._Tsub_mean.tif'};
        
%%

%% Set options
options.maxVal = 2^13;              % Maximum pixel intensity. Set anything above to 0

% Size of image patches
options.patchSize = [64 64];        % Size of image patches
options.stride = 56;                % Offset of image patches

% Initial ROI detection parameters
options.DETREND_FRAMES = 45;        % Number of frames to compute baseline F0 over
options.thr = 3;                    % Threshold for active pixels 
options.filtSize = 3;               % Median filter size for initial ROI detection
options.overlapThr = 0.5;           % Spatial overlap merge threshold 
options.sizeRange = [30 2000];      % Allowable size range of valid ROIs
options.eta = 0.01;                 % Temporal regularization weight
options.beta = 0.5;                 % Spatial regularization weight
options.med_opt = false;

% ROI cleanup parameters
options.thr_method = 'quant';       % Method of thresholding ROIs: 'max' or 'quant'
options.quantileThr = 0.9;          % Quantile threshold for thresholding using 'quant'
options.maxthr = 0.2;               % Max threshold for thresholding using 'max'
options.final_C = true;             % Whether or not to recompute traces C after ROI cleanup

% ROI validity & merging parameters
options.minSkew = 0;                % Minimum skew of temporal trace of valid ROIs
options.shapeThr = 0.5;             % Correlation threshold for ROIs
options.temporalCorrThr = 0.9;     % Temporal correlation merge threshold 


tempFolder = '.';
for i_file = length(files)
    thisFile = files{i_file};
    
    [cROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(thisFile, options);
    
    folder = fileparts(thisFile);
    outFolder = folder;
    if(~exist(outFolder,'dir'))
        mkdir(outFolder);
    end
    fprintf('Saving %d ROIs...',size(cROIs,2));
    save(fullfile(outFolder,'DNMF_Out.mat'), 'cROIs', 'Cs', 'coherence', 'skew', 'sz', 'options', 'tElapsed', '-v7.3');
        
    if(~isempty(tempFolder))
        outFolder_original = outFolder;
        outFolder = tempFolder;
    end
    outputFolder = fullfile(outFolder,'ROI_Set\');
    if(~exist(outputFolder,'dir'))
        mkdir(outputFolder);
    end
    [~,max_idx] = max(Cs,[],2);
    for i_roi = 1:size(cROIs,2)   
        temp = reshape(full(cROIs(:,i_roi)),[512 512]);
        [a,b] = find(imdilate(temp>0,ones(3)));
        c = boundary(a,b,0.95);

        writeImageJROI_3([a(c) b(c)], 4, max_idx(i_roi), sprintf('r%04d',i_roi), outputFolder);
    end

    zip(strrep(outputFolder,'ROI_Set\','RoiSet_Auto'),'*',outputFolder);
    rmdir(outputFolder, 's');
    if(~isempty(tempFolder))
        movefile(fullfile(outFolder,'RoiSet_Auto.zip'),fullfile(outFolder_original,'RoiSet_Auto.zip'));        
    end
end