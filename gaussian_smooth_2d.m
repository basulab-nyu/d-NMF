function [newDensity] = gaussian_smooth_2d(oldDensity,sig,D)
    %[newDensity] = gaussian_smooth_2d(oldDensity,sig,D)    
    if(nargin<2)
        sig = 1;
    end
    
    if(nargin<3)
        D = ceil(3*sig*[1 1]);
    end
    
    newDensity = oldDensity;
    if(sig>0)
        g = fspecial('gaussian',D,sig);        
        newDensity = convolve2(newDensity,g,'replicate');                    
    end
end