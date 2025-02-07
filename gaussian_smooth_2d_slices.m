function [sm_data] = gaussian_smooth_2d_slices(data,sig,D)

    sm_data = NaN(size(data));

    if(nargin<2)
        sig = 1;
    end
    
    if(nargin<3)
        D = ceil(3*sig*[1 1]);
    end   

    for i=1:size(data,3)
        sm_data(:,:,i) = gaussian_smooth_2d(data(:,:,i),sig,D);        
    end

end