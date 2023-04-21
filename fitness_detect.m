function [START_STOP, PEAK] = fitness_detect(traces_C, rois_A, video_Y)
    % [START_STOP, PEAK] = fitness_detect(traces_C, rois_A, video_Y)
    
    THR_DF = 1.95;
    THR_FIT = 0.28;
    DETREND_FRAMES = 45;
    DT = 20/30;

    video_Y = j_detrend2b(video_Y,DETREND_FRAMES,3);
    video_Y(isnan(video_Y) | isinf(video_Y)) = 0;
    zy = zscore(video_Y,[],3);    
    
    
    
    THIS1 = traces_C';
    temp = imerode(THIS1,ones([DETREND_FRAMES,1]));
    THIS1 = zscore(THIS1-temp)';

    zCCC = corrTrace(zy, rois_A, 3);
    METHOD = @(x,y) fp_detect(x, DT, y);
    
    START_STOP = cell(1,size(rois_A,3));
    PEAK = cell(1,size(rois_A,3));
    
    for i_roi = 1:size(rois_A,3)
        zcs2 = THIS1(i_roi,:);
        [startStop, peakIndex, amp, width] = METHOD(zcs2,THR_DF);
        THIS2 = zCCC(i_roi,:);
        screenVals = THIS2(peakIndex);
        valid = screenVals>THR_FIT;
        START_STOP{i_roi} = startStop(valid,:);
        PEAK{i_roi} = peakIndex(valid);
    end     

end