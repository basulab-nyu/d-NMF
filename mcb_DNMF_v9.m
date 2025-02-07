function [cROIs, Cs, cROIs_BG, Cs_BG, coherence, skew, sz, patchID, tElapsed, dimensions] = mcb_DNMF_v9(paths_to_videos, options)
    % [cROIs, Cs, cROIs_BG, Cs_BG, coherence, skew, sz, tElapsed, dimensions] = mcb_DNMF_v9(path_to_video, options)
    %
    % Jason Moore, 2025
    
    
    
    % Load in Video
    fileName = paths_to_videos{1};
    Y = bigread2(fileName,1,1);    
    
    tStart = tic;
    [szA, szB, ~] = size(Y);
    
    if(nargin<2)
        options = defaultOptions_mcbDNMF();
    end
    
    minSkew = options.minSkew;
    sizeRange = options.sizeRange;
    options.szA = szA;
    options.szB = szB;
    
    % Detect ROIs in patches
    [cROIs0, Cs0, cROIs_BG, Cs_BG, coherence0, skew0, sz0, patchID0] = DNMF_General9(paths_to_videos, options);
    
    % Light screen for skewness and ROI size
    valid = skew0>minSkew & isbetween(sz0,sizeRange(1),sizeRange(2));
    cROIs = cROIs0(:,valid);
    Cs = Cs0(valid,:);
    coherence = coherence0(valid);
    skew = skew0(valid);
    sz = sz0(valid);
    patchID = patchID0(valid);
    dimensions = [szA szB];
    
    %% No merging. Leave that to the user later
    tElapsed = toc(tStart);
    fprintf('Done!\n');
end