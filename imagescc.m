function [] = imagescc(c)
    % [] = imagescc(c)
    % Reshapes a column vector into a square array and then calls imagesc
    % For rapid viewing of ROIs
    
    d = size(c,1);
    imagesc(reshape(c,[sqrt(d) sqrt(d)]));
end