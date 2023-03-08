function [startStopIndices] = j_event_detect(Y, onThr, offThr)
    % [startStopIndices] = j_event_detect(Y, onThr, offThr)

    Y = zscore(Y);
    startIndices = jzeroCrossing(Y-onThr,1);
    if(isempty(startIndices))
        startStopIndices = NaN(0,2);
        return;
    end
    stopIndices = jzeroCrossing(Y-offThr,-1);
    if(isempty(stopIndices) && ~isempty(startIndices))
        stopIndices = length(Y);
    end
    
    if(startIndices(end)>stopIndices(end))
        stopIndices = [stopIndices length(Y)];
    end
        
    id = nearestpoint(startIndices,stopIndices,'next');
    stopIndices = stopIndices(id);
    
    [b,i,j] = unique(stopIndices);

    startIndices = startIndices(i);
    stopIndices = stopIndices(i);
    
    startStopIndices = [startIndices(:) stopIndices(:)];       
end