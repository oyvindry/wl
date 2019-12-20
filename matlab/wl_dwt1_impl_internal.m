function x=wl_dwt1_impl_internal(x, f, prefilter, offsets, m, bd_mode, data_layout)
    % Compute a 1D DWT using a precomputed kernel. The kernel may be the default library kernel obtained by calling find_kernel, 
    % or may be used-defined.
    %
    % x:         Matrix whose DWT will be computed along the first dimension(s). 
    % f:         kernel function     
    % prefilter: function which computes prefiltering. The default is no prefiltering.
    % offsets:   offsets at the beginning and the end as used by boundary wavelets. 
    % m:              Number of resolutions. 
    % bd_mode:        Boundary extension mode. Possible modes are. 
    %                 'per'    - Periodic extension 
    %                 'symm'   - Symmetric extension 
    %                 'none'   - Take no extra action at the boundaries (i.e., the boundaries are zero-padded)
    %                 'bd'     - Preserve vanishing moments at the boundaries
    % data_layout:    How data should be assembled. Possible modes are:
    %                 'resolution': Lowest resolution first
    %                 'time': Sort according to time

    inds = 1:size(x,1);
    x = prefilter(x, 1);
    for res=0:(m - 1)
        x(inds, :) = f(x(inds, :), bd_mode);
        inds = inds((offsets(1,1)+1):2:(end-offsets(1,2)));
    end
    % x(inds, :) = prefilter(x(inds, :),0);
    x = reorganize_coeffs_forward(x, m, offsets, data_layout);
end

function y=reorganize_coeffs_forward(x, m, offsets, data_layout)
    y = x;
    if strcmpi(data_layout, 'resolution')
        N = size(x,1);
        inds = 1:N;
        endy = N;
        for res=1:m
            xindices = [inds(1:offsets(1,1)) inds((offsets(1,1) + 2):2:(end-offsets(1,2))) inds((end-offsets(1,2)+1):end)]; % psi-indices
            y((endy-length(xindices)+1):endy,:) = x(xindices,:);
            endy = endy-length(xindices);
            inds = inds((offsets(1,1)+1):2:(end-offsets(1,2))); 
        end
        y(1:endy, :) = x(inds, :);
    end
end
