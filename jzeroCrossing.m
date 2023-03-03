function [points] = jzeroCrossing(signal,direction)
    % [points] = jzeroCrossing(signal,direction)
    
    if(nargin<2)
        direction = 0;
    end

    signal(signal==0) = 1e-10;
    
    %signal = reshape(signal,length(signal),[]);
    side = sign(signal);
    ddd = [0 reshape(diff(side),1,[])];
    if(direction==0)
        points = find(abs(ddd)>0);
    elseif(direction==1)
        points = find(ddd>0);
    elseif(direction==-1)
        points = find(ddd<0);
    else
        fprintf('Invalid direction, returning empty\n');
        points = [];
    end
end