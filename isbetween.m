function [out] = isbetween(numbers, low, high)
    % [out] = isbetween(numbers, low, high)
    if(isempty(low))
        out = (zeros(size(numbers))==1)';        
    else
        out = numbers>=low & numbers<=high;
    end
end