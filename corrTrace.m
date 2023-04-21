function [corTrace] = corrTrace(Y, ROIs, padding)
    % [corTrace] = corrTrace(Y, ROIs, padding)
    
    if(nargin<3)
        padding = 3;
    end
    corTrace = NaN(size(ROIs,3),size(Y,3));
    for i_roi = 1:size(ROIs,3)
        this = ROIs(:,:,i_roi);
        [b,a] = find(this>0);
        bb = [min(b) max(b)]+padding*[-1 1];
        aa = [min(a) max(a)]+padding*[-1 1];
        bb(bb<1) = 1;
        bb(bb>size(Y,1)) = size(Y,1);
        aa(aa<1) = 1;
        aa(aa>size(Y,2)) = size(Y,2);
        cThis = reshape(this(bb(1):bb(2),aa(1):aa(2)),[],1);
        cY = reshape(Y(bb(1):bb(2),aa(1):aa(2),:),[],size(Y,3));
        C = zscore(cThis)'*zscore(cY)/(length(cThis)-1);
        corTrace(i_roi,:) = C;
    end


end