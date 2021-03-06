function x=wl_idwt3_impl_internal(x, fx, fy, fz, prefilterx, prefiltery, prefilterz, offsets, varargin)
    % Compute a 3D IDWT using precomputed kernels. The kernels may be the default library kernels obtained by calling find_kernel, 
    % or may be used-defined.
    %
    % x:         Matrix whose IDWT3 will be computed along the first dimensions. 
    % fx, fy, fz: kernel functions     
    % prefilterx, prefiltery, prefilterz: functions which compute prefiltering. The default is no prefiltering.
    % offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    %
    % This function also accepts a number of named, optional parameters. These are parsed by the function wl_setopts(). 
    % The documentation of this function also contains the full documentation of these optional parameters.
    
    opts = wl_setopts('wave_name', 'unknown', 'dims', 3, varargin{:});
    
    [x, resstart, resend] = reorganize_coeffs3_reverse(x, opts.m, offsets, opts.data_layout);
    
    % postconditioning
    indsx = resstart(1,opts.m+1):2^opts.m:resend(1,opts.m+1); indsy = resstart(2,opts.m+1):2^opts.m:resend(2,opts.m+1); indsz = resstart(3,opts.m+1):2^opts.m:resend(3,opts.m+1);
    % x(indsx, indsy, indsz, :) = tensor3_impl(x(indsx, indsy, indsz, :), @(x,bdm) prefilterx(x, 1), @(x,bdm) prefiltery(x, 1), @(x,bdm) prefilterz(x, 1), opts.bd_mode);

    for res = (opts.m - 1):(-1):0
        indsx = resstart(1,res+1):2^res:resend(1,res+1); 
        indsy = resstart(2,res+1):2^res:resend(2,res+1);
        indsz = resstart(3,res+1):2^res:resend(3,res+1);
        x(indsx, indsy, indsz, :) = tensor3_impl(x(indsx, indsy, indsz, :), fx, fy, fz, opts.bd_mode);
    end
    
    % preconditioning
    indsx = resstart(1,1):resend(1,1); indsy = resstart(2,1):resend(2,1); indsz = resstart(3,1):resend(3,1);
    x(indsx, indsy, indsz, :) = tensor3_impl(x(indsx, indsy, indsz, :), @(x,bdm) prefilterx(x, 0), @(x,bdm) prefiltery(x, 0), @(x,bdm) prefilterz(x, 0), opts.bd_mode);
end

function [sig_out, resstart, resend]=reorganize_coeffs3_reverse(sig_in, m, offsets, data_layout)
    indsx = 1:size(sig_in,1); indsy = 1:size(sig_in,2); indsz = 1:size(sig_in,3);
    sig_out = sig_in;
    resstart = [1:(m+1); 1:(m+1); 1:(m+1)]; resend = [1:(m+1); 1:(m+1); 1:(m+1)];
    resstart(1,1) = indsx(1); resend(1,1) = indsx(end);
    resstart(2,1) = indsy(1); resend(2,1) = indsy(end);
    resstart(3,1) = indsz(1); resend(3,1) = indsz(end);
    if strcmpi(data_layout, 'time')
        for res=1:m
            indsx = indsx((offsets(1,1)+1):2:(end-offsets(1,2)));
            indsy = indsx((offsets(2,1)+1):2:(end-offsets(2,2)));
            indsz = indsz((offsets(3,1)+1):2:(end-offsets(3,2)));
            resstart(1,res+1) = indsx(1); resend(1,res+1) = indsx(end);
            resstart(2,res+1) = indsy(1); resend(2,res+1) = indsy(end);
            resstart(3,res+1) = indsz(1); resend(3,res+1) = indsz(end);
        end
    end
    if strcmpi(data_layout, 'resolution')
        endx = size(sig_in,1); endy = size(sig_in,2); endz = size(sig_in,3);
        for res=1:m
            psiinds_x = [indsx(1:offsets(1,1)) indsx((offsets(1,1) + 2):2:(end-offsets(1,2))) indsx((end-offsets(1,2)+1):end)]; % psi-indices
            psiinds_y = [indsy(1:offsets(2,1)) indsy((offsets(2,1) + 2):2:(end-offsets(2,2))) indsy((end-offsets(2,2)+1):end)];
            psiinds_z = [indsz(1:offsets(3,1)) indsz((offsets(3,1) + 2):2:(end-offsets(3,2))) indsz((end-offsets(3,2)+1):end)];
            phiinds_x = indsx((offsets(1,1) + 1):2:(end-offsets(1,2)));
            phiinds_y = indsx((offsets(2,1) + 1):2:(end-offsets(2,2)));
            
            resstart(1,res+1) = indsx(offsets(1,1)+1); resend(1,res+1)   = indsx(end-offsets(1,2));
            resstart(2,res+1) = indsy(offsets(2,1)+1); resend(2,res+1)   = indsy(end-offsets(2,2));
            resstart(3,res+1) = indsz(offsets(3,1)+1); resend(3,res+1)   = indsz(end-offsets(3,2));
            
            sig_out(psiinds_x, indsy, indsz, :) = sig_in((endx-length(psiinds_x)+1):endx, 1:endy, 1:endz, :);
            sig_out(phiinds_x, psiinds_y, indsz, :) = sig_in(1:(endx-length(psiinds_x)), (endy-length(psiinds_y)+1):endy, 1:endz, :);
            sig_out(phiinds_x, phiinds_y, psiinds_z, :) = sig_in(1:(endx-length(psiinds_x)), 1:(endy-length(psiinds_y)), (endz-length(psiinds_z)+1):endx, :);
            
            endx = endx - length(psiinds_x); endy = endy - length(psiinds_y); endz = endz - length(psiinds_z);
            indsx = indsx((offsets(1,1)+1):2:(end-offsets(1,2))); 
            indsy = indsy((offsets(2,1)+1):2:(end-offsets(2,2)));
            indsz = indsz((offsets(3,1)+1):2:(end-offsets(3,2)));
        end
        sig_out(indsx, indsy, indsz, :) = sig_in(1:endx, 1:endy, 1:endz, :);
    end
end
