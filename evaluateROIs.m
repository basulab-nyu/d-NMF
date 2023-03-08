function [coherence, skew, roiSize] = evaluateROIs(A, C, patchSize)
    % [coherence, skew, roiSize] = evaluateROIs(A, C, patchSize)
    % Computes coherence, skew, and size of an ROI
    % Coherence is computed on all pixels in the FOV, so small ROIs will be
    % biased towards high values due to most pixels being 0
    
    A2 = reshape(full(A), patchSize(1), patchSize(2), size(A,2));
    h = [1 1 1; 1 0 1; 1 1 1];
    coherence = zeros(size(C,1),1);
    for i_a = 1:size(A2,3)
        this = squeeze(A2(:,:,i_a));
        c = corrcoef(this, conv2(this,h,'same'));
        coherence(i_a) = c(1,2);
    end
    skew = skewness(full(C),[],2);

    roiSize = full(sum(A>0))';
end