function [cROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(path_to_video, options)
    % [cROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(path_to_video, options)
    
    % Load in Video
    fileName = path_to_video;
    Y = bigread2(fileName,1);
    
    % Set saturated pixels to 0
    Y(Y>options.maxVal) = 0;
    
    tStart = tic;
    [szA, szB, ~] = size(Y);
    
    if(nargin<2)
        options = defaultOptions_mcbDNMF();
    end
    
    temporalCorrThr = options.temporalCorrThr;
    minSkew = options.minSkew;
    sizeRange = options.sizeRange;
        
    
    % Detect ROIs in patches
    [cROIs0, Cs0, coherence0, skew0, sz0] = DNMF_General3(Y, options);
    
    valid = skew0>minSkew & isbetween(sz0,sizeRange(1),sizeRange(2));
    cROIs1 = cROIs0(:,valid);
    Cs1 = Cs0(valid,:);
    
    
    
    %% Merge ROIs
    fprintf('Merging ROIs...');
    
    % De-trend and z-score to compute temporal correlation
    Cs1b = (j_detrend2b(Cs1,options.DETREND_FRAMES,2,false));
    zCs = zscore(Cs1b,[],2);
    
    % Compute spatial overlap
    cROIs_binary = double(cROIs1>0);
    spatialOverlap_0 = cROIs_binary'*cROIs_binary;
    sizes = sum(cROIs_binary);
    [sz1, sz2] = meshgrid(sizes);
    temp = zeros(size(sz1,1),size(sz1,2),2);
    temp(:,:,1) = sz1;
    temp(:,:,2) = sz2;
    minSize = min(temp,[],3);

    spatialOverlap = spatialOverlap_0./minSize;
    temporalCorr = zCs*zCs'/(size(zCs,2)-1);
    
    % Find connected components of the graph
    linked = (spatialOverlap>0) & (temporalCorr>temporalCorrThr);    
    [~,connComp] = graphconncomp(sparse(linked));

    % Merge ROIs that are part of connected components
    cROIs = sparse(szA*szB, max(connComp));
    Cs = zeros(max(connComp), size(Cs1,2));
    coherence = NaN(max(connComp),1);
    skew = NaN(max(connComp),1);
    sz = NaN(max(connComp),1);
    for i_component=1:max(connComp)
        these = connComp==i_component; 
        if(sum(these)==1)                   % Connected components with a single member
            cROIs(:,i_component) = cROIs1(:,these);
            Cs(i_component,:) = Cs1(these,:);
        else                                % Connected components with multiple members
            rois = cROIs1(:,these);
            traces = Cs1(these,:);
            maxTraces = max(traces,[],2);
            rois2 = bsxfun(@times,rois,maxTraces');      % Scale ROIs by the maximum value of their temporal trace
            temp = max(rois2,[],2);
            n = sqrt(sum(temp.^2));
            temp = temp/n;                  % New ROI is the maximum of member ROIs, normalized to have length 1
            cROIs(:,i_component) = temp;
            
            indices = find(these);
            validPixels = temp>0;
            values = zeros(sum(validPixels),size(Cs,2));
            for i_part = 1:length(indices)
                values = max(values, rois(validPixels,i_part)*traces(i_part,:));                
            end
            traces = temp(validPixels)\values;        % Solve for temporal component that best represents the new ROI
            
            
            Cs(i_component,:) = traces;
        end
        [ch, sk, z] = evaluateROIs(cROIs(:,i_component), Cs(i_component,:),[szA szB]);
        coherence(i_component) = ch;
        skew(i_component) = sk;
        sz(i_component) = z;
    end
   
    tElapsed = toc(tStart);
    fprintf('Done!\n');
end