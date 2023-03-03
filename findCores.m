function [ROI] = findCores(V, thr, filtMask, minSz)

    if(nargin<4)
        minSz = 15;
    end

    if(nargin<3)
        filtMask = [3 3 3];
    end
    if(nargin<2)
        thr = 2;
    end

    height = size(V,1);
    width = size(V,2);

    CC = bwconncomp(medfilt3(bsxfun(@gt,V,thr),filtMask));
    sz = cellfun(@(x) length(x), CC.PixelIdxList);
    segs = CC.PixelIdxList(sz>2*minSz);

    ROIs = zeros(height,width,length(segs));
    for jj=1:length(segs)
        temp = mod(segs{jj}-1,height*width)+1;
        temp2 = zeros(height,width);
        temp2(temp) = 1;
        ROIs(:,:,jj) = temp2;    
    end
    sz = squeeze(sum(sum(ROIs,1),2));
    valid = sz>minSz;
    ROI = ROIs(:,:,valid);
    
end