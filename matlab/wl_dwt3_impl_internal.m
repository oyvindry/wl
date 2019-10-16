function x = wl_dwt3_impl_internal(x, fx, fy, fz, prefilterx, prefiltery, prefilterz, offsets, varargin)
    % Compute a 3D DWT using precomputed kernels. The kernels may be the default library kernels obtained by calling find_kernel, 
    % or may be used-defined.
    %
    % x:         Matrix whose DWT3 will be computed along the first dimensions. 
    % fx, fy, fz: kernel functions     
    % prefilterx, prefiltery, prefilterz: functions which compute prefiltering. The default is no prefiltering.
    % offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    %
    % This function also accepts a number of named, optional parameters. These are parsed by the function wl_setopts(). 
    % The documentation of this function also contains the full documentation of these optional parameters.
    
    opts = wl_setopts('wave_name', 'unknown', 'dims', 3, varargin{:});
    
    indsx = 1:size(x,1); indsy = 1:size(x,2); indsz = 1:size(x,3);

    % preconditioning   
    x(indsx, indsy, indsz, :) = tensor3_impl(x(indsx, indsy, indsz, :), @(x,bdm) prefilterx(x, 1), @(x,bdm) prefiltery(x, 1), @(x,bdm) prefilterz(x, 1), opts.bd_mode);
         
    for res = 0:(opts.m - 1)
        x(indsx, indsy, indsz, :) = tensor3_impl(x(indsx, indsy, indsz, :), fx, fy, fz, opts.bd_mode);  
        indsx = indsx((offsets(1,1)+1):2:(end-offsets(1,2)));   
        indsy = indsy((offsets(2,1)+1):2:(end-offsets(2,2)));
        indsz = indsz((offsets(3,1)+1):2:(end-offsets(3,2)));
    end
    
    % postconditioning
    x(indsx, indsy, indsz,:) = tensor3_impl(x(indsx, indsy, indsz, :), @(x,bdm) prefilterx(x, 0), @(x,bdm) prefiltery(x, 0), @(x,bdm) prefilterz(x, 0), opts.bd_mode);  
    
    x = reorganize_coeffs3_forward(x, opts.m, offsets, opts.data_layout);   
end

function sig_out=reorganize_coeffs3_forward(sig_in, m, offsets, data_layout)
    sig_out = sig_in;
    if strcmpi(data_layout, 'resolution')
        indsx = 1:size(sig_in,1); indsy = 1:size(sig_in,2); indsz = 1:size(sig_in,3);
        endx = size(sig_in,1); endy = size(sig_in,2); endz = size(sig_in,3); 
        for res=1:m
            psiinds_x = [indsx(1:offsets(1,1)) indsx((offsets(1,1) + 2):2:(end-offsets(1,2))) indsx((end-offsets(1,2)+1):end)]; % psi-indices
            psiinds_y = [indsy(1:offsets(2,1)) indsy((offsets(2,1) + 2):2:(end-offsets(2,2))) indsy((end-offsets(2,2)+1):end)];
            psiinds_z = [indsz(1:offsets(3,1)) indsz((offsets(3,1) + 2):2:(end-offsets(3,2))) indsz((end-offsets(3,2)+1):end)];
            phiinds_x = indsx((offsets(1,1) + 1):2:(end-offsets(1,2)));
            phiinds_y = indsy((offsets(2,1) + 1):2:(end-offsets(2,2)));
            
            sig_out( (endx-length(psiinds_x)+1):endx, 1:endy, 1:endz, :) = sig_in(psiinds_x, indsy, indsz, :);
            sig_out( 1:(endx-length(psiinds_x)), (endy-length(psiinds_y)+1):endy, 1:endz, :) = sig_in(phiinds_x, psiinds_y, indsz, :);
            sig_out( 1:(endx-length(psiinds_x)), 1:(endy-length(psiinds_y)), (endz-length(psiinds_z)+1):endx, :) = sig_in(phiinds_x, phiinds_y, psiinds_z, :);
            
            endx = endx - length(psiinds_x); endy = endy - length(psiinds_y); endz = endz - length(psiinds_z);
            indsx = indsx((offsets(1,1)+1):2:(end-offsets(1,2))); 
            indsy = indsy((offsets(2,1)+1):2:(end-offsets(2,2)));
            indsz = indsz((offsets(3,1)+1):2:(end-offsets(3,2)));
        end
        sig_out(1:endx, 1:endy, 1:endz, :) = sig_in(indsx, indsy, indsz, :);
    end
end