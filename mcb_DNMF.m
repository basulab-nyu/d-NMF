function [cROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(path_to_video, options)
    % [ROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(path_to_video, options)
    fileName = path_to_video;
    Y = bigread2(fileName,1);
%     Y = double(Y);
    Y(Y>options.maxVal) = 0;
    
    tStart = tic;
    [szA, szB, nFrames] = size(Y);
    
    if(nargin<2)
        options = defaultOptions_mcbDNMF();
    end
    
    thr = options.thr;
    patchSize = options.patchSize;
    stride = options.stride;
    overlapThr = options.overlapThr;
    temporalCorrThr = options.temporalCorrThr;
    minSkew = options.minSkew;
    sizeRange = options.sizeRange;
        
%     [cROIs0, Cs0, coherence0, skew0, sz0] = DNMF_General3(Y, thr, patchSize, stride, overlapThr, sizeRange);
    [cROIs0, Cs0, coherence0, skew0, sz0] = DNMF_General3(Y, options);
    
    fprintf('Merging ROIs...');
    valid = skew0>minSkew & isbetween(sz0,sizeRange(1),sizeRange(2));
    cROIs1 = cROIs0(:,valid);
    Cs1 = Cs0(valid,:);
    coherence1 = coherence0(valid);
    skew1 = skew0(valid);
    sz1 = sz0(valid);
    
%     ROIs1 = reshape(full(cROIs1),[szA, szB, size(cROIs1,2)]);    

    
    %%
    Cs1b = (j_detrend2b(Cs1,options.DETREND_FRAMES,2,false));
    zCs = zscore(Cs1b,[],2);
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
    
    linked = (spatialOverlap>0) & (temporalCorr>temporalCorrThr);    
    [SSS,CCC] = graphconncomp(sparse(linked));

    
    cROIs = sparse(szA*szB, max(CCC));
    Cs = zeros(max(CCC), size(Cs1,2));
    coherence = NaN(max(CCC),1);
    skew = NaN(max(CCC),1);
    sz = NaN(max(CCC),1);
    for ii=1:max(CCC)
        these = CCC==ii; 
        if(sum(these)>1)
            r = cROIs1(:,these);
%             r = reshape(ROIs1(:,:,these),512*512,[]);
            c = Cs1(these,:);
            mc = max(c,[],2);
            r2 = bsxfun(@times,r,mc');
            temp = max(r2,[],2);
            n = sqrt(sum(temp.^2));
            temp = temp/n;
            q = find(these);
            validPixels = temp>0;
            v = zeros(sum(validPixels),size(Cs,2));
            for i_part = 1:length(q)
                v = max(v, r(validPixels,i_part)*c(i_part,:));                
            end
            c = temp(validPixels)\v;
            
            cROIs(:,ii) = temp;
%             ROIs(:,:,ii) = reshape(temp,[szA szB]);
            Cs(ii,:) = c;
        else
            cROIs(:,ii) = cROIs1(:,these);
%             ROIs(:,:,ii) = ROIs1(:,:,these);
            Cs(ii,:) = Cs1(these,:);
        end
        [ch, sk, z] = evaluateROIs(cROIs(:,ii), Cs(ii,:),[szA szB]);
        coherence(ii) = ch;
        skew(ii) = sk;
        sz(ii) = z;
    end
   
    tElapsed = toc(tStart);
    fprintf('Done!\n');
end