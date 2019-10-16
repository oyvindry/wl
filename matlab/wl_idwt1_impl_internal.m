function x=wl_idwt1_impl_internal(x, f, prefilter, offsets, varargin)
    % Compute a 1D IDWT using a precomputed kernel. The kernel may be the default library kernel obtained by calling find_kernel, 
    % or may be used-defined.
    %
    % x:         Matrix whose IDWT will be computed along the first dimension(s). 
    % f:         kernel function     
    % prefilter: function which computes prefiltering. The default is no prefiltering.
    % offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    %
    % This function also accepts a number of named, optional parameters. These are parsed by the function wl_setopts(). 
    % The documentation of this function also contains the full documentation of these optional parameters.
    
    opts = wl_setopts('wave_name', 'unknown', 'dims', 1, varargin{:});
    
    [x, resstart, resend] = reorganize_coeffs_reverse(x, opts.m, offsets, opts.data_layout);
    x(resstart(opts.m+1):2^opts.m:resend(opts.m+1),:) = prefilter(x(resstart(opts.m+1):2^opts.m:resend(opts.m+1),:), 1);
    for res = (opts.m - 1):(-1):0
        inds = resstart(res+1):2^res:resend(res+1);
        x(inds, :) = f(x(inds, :), opts.bd_mode);
    end
    x = prefilter(x, 0);
end

function [y, resstart, resend]=reorganize_coeffs_reverse(x, m, offsets, data_layout)
    N = size(x,1);
    inds = 1:N;
    y = x;
    resstart = 1:(m+1); resend = 1:(m+1);
    resstart(1) = inds(1);
    resend(1) = inds(end);
    if strcmpi(data_layout, 'time')
        for res=1:m
            inds = inds((offsets(1,1)+1):2:(end-offsets(1,2)));
            resstart(res+1) = inds(1);
            resend(res+1) = inds(end);
        end
    end
    if strcmpi(data_layout, 'resolution')
        endy = N;
        for res=1:m
            xindices = [inds(1:offsets(1,1)) inds((offsets(1,1) + 2):2:(end-offsets(1,2))) inds((end-offsets(1,2)+1):end)]; % psi-indices
            resstart(res+1) = inds(offsets(1,1)+1);
            resend(res+1)   = inds(end-offsets(1,2));
            y(xindices,:) = x((endy-length(xindices)+1):endy,:);
            endy = endy-length(xindices);
            inds = inds((offsets(1,1)+1):2:(end-offsets(1,2)));
        end
        y(inds,:) = x(1:endy,:);
    end
end