function [] = imagescc(c, dimensions)
    % [] = imagescc(c, dimensions)
    % Reshapes a column vector into a square array and then calls imagesc
    % For rapid viewing of ROIs
    
    if(nargin<2)
        d = size(c,1);
        dimensions = [sqrt(d) sqrt(d)];
    end
    imagesc(reshape(c,dimensions));
end