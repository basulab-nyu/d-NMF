function [Y_out, Y0] = j_detrend2b(Y_in, L, dim, norm_opt, med_opt)
    % [Y_out, Y0] = j_detrend2b(Y_in, L, dim, norm_opt, med_opt)
    % Input
    % L: Number of frames over which to compute moving baseline Y0
    % dim: Dimension across which to compute Y0
    % norm_opt: Option to compute output as (Y-Y0)/Y0 (true) or (Y-Y0)
    % (false)
    % med_opt: Option to compute Y0 as rolling median (true) or minimum
    % (false)
    % Output
    % Y_out: Detrended version of Y_in
    % Y0: Moving baseline
    %
    % Jason Moore, 2025
    if(nargin<5)
        med_opt = false;
    end
    if(nargin<4)
        norm_opt = true;
    end
    if(nargin<3)
        dim = 3;
    end
    if(L<=1)
        Y_out = Y_in;
        return;
    end
    
    % Compute rolling Y0
    d = ones(1, dim);
    d(dim) = L;
    if(med_opt)     % Median
        Y0 = movmedian(Y_in, L, dim);
    else            % Minimum
        Y0 = imerode(Y_in, ones(d));
    %     Y0 = movmin(Y_in, L, dim);
    end
    if(norm_opt)
        Y_out = (Y_in - Y0)./Y0;    % Subtract and divide if normalizing
    else
        Y_out = (Y_in - Y0);        % Just subtract if only detrending
    end
end