function [eventStartStop, eventPeakLocation, eventAmp, eventWidth] = fp_detect(Y, dt, minHeight)
    % [eventStartStop, eventPeakLocation, eventAmp, eventWidth] = fp_detect(Y, dt, minHeight)
    
    Y = zscore(Y);
    startStopPositive = j_event_detect(Y, minHeight, 0.5);
    startStopNegative = j_event_detect(-Y, minHeight, 0.5);
    startStopPositive(:,2) = max([startStopPositive(:,1) startStopPositive(:,2)-1],[],2);
    [ampsPositive,peakLocPositive] = arrayfun(@(x) max(Y(startStopPositive(x,1):startStopPositive(x,2))), 1:size(startStopPositive,1));
    [ampsNegative,peakLocNegative] = arrayfun(@(x) max(-Y(startStopNegative(x,1):startStopNegative(x,2))), 1:size(startStopNegative,1));
    
    if(isempty(ampsPositive))
        eventStartStop = NaN(0,2);
        eventPeakLocation = NaN(0,0);
        eventAmp = NaN(0,0);
        eventWidth = NaN(0,0);
        return;
    end
    widthsPositive = diff(startStopPositive,[],2)*dt;
    widthsNegative = diff(startStopNegative,[],2)*dt;

    edgesAmp = 0:0.5:max([ampsPositive(:);ampsNegative(:)])+0.5;
    edgesWidth = 0:0.25:max([widthsPositive(:);widthsNegative(:)])+0.25;
    
    [countPositive,~,~,binPositive] = histcn([ampsPositive(:) widthsPositive(:)],edgesAmp, edgesWidth);
    idxPositive = sub2ind(size(countPositive),binPositive(:,1),binPositive(:,2));    
    [countNegative,~,~,binNegative] = histcn([ampsNegative(:) widthsNegative(:)],edgesAmp, edgesWidth);
    idxNegative = sub2ind(size(countNegative),binNegative(:,1),binNegative(:,2));    
    
    fpRate = countNegative./countPositive;
    valid = fpRate<0.05;
    validIndices = find(valid);
    
    validEvents = ismember(idxPositive, validIndices);
   
    eventStartStop = startStopPositive(validEvents,:);
    eventAmp = ampsPositive(validEvents);
    eventWidth = widthsPositive(validEvents);
    eventPeakLocation = startStopPositive(validEvents,1)+peakLocPositive(validEvents)'-1;
end