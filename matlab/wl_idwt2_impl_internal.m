function x=idwt2_impl_internal(x, fx, fy, prefilterx, prefiltery, offsets,  m, bd_mode, data_layout)
    % Compute a 2D IDWT using precomputed kernels. The kernels may be the default library kernels obtained by calling find_kernel, 
    % or may be used-defined.
    %
    % x:         Matrix whose IDWT2 will be computed along the first dimensions. 
    % fx, fy:    kernel functions     
    % prefilterx, prefiltery: functions which compute prefiltering. The default is no prefiltering.
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
    
    [x, resstart, resend] = reorganize_coeffs2_reverse(x, m, offsets, data_layout);
    % postconditioning
    indsx = resstart(1,m+1):2^m:resend(1,m+1); indsy = resstart(2,m+1):2^m:resend(2,m+1);
    % x(indsx, indsy, :) = tensor2_impl(x(indsx, indsy, :), @(x,bdm) prefilterx(x, 1), @(x,bdm) prefiltery(x, 1), bd_mode);

    for res = (m - 1):(-1):0
        indsx = resstart(1,res+1):2^res:resend(1,res+1); 
        indsy = resstart(2,res+1):2^res:resend(2,res+1);
        x(indsx, indsy, :) = tensor2_impl(x(indsx, indsy, :), fx, fy, bd_mode);
    end
    
    % preconditioning
    indsx = resstart(1,1):resend(1,1); indsy = resstart(2,1):resend(2,1);
    x(indsx, indsy, :) = tensor2_impl(x(indsx, indsy, :), @(x,bdm) prefilterx(x, 0), @(x,bdm) prefiltery(x, 0), bd_mode);
end

function [sig_out, resstart, resend]=reorganize_coeffs2_reverse(sig_in, m, offsets, data_layout)
    indsx = 1:size(sig_in,1); indsy = 1:size(sig_in,2);
    sig_out = sig_in;
    resstart = [1:(m+1); 1:(m+1)]; resend = [1:(m+1); 1:(m+1)];
    resstart(1,1) = indsx(1); resend(1,1) = indsx(end);
    resstart(2,1) = indsy(1); resend(2,1) = indsy(end);
    if strcmpi(data_layout, 'time')
        for res=1:m
            indsx = indsx((offsets(1,1)+1):2:(end-offsets(1,2)));
            indsy = indsy((offsets(2,1)+1):2:(end-offsets(2,2)));
            resstart(1,res+1) = indsx(1); resend(1,res+1) = indsx(end);
            resstart(2,res+1) = indsy(1); resend(2,res+1) = indsy(end);
        end
    end
    if strcmpi(data_layout, 'resolution')
        endx = size(sig_in,1); endy = size(sig_in,2);
        for res=1:m
            psiinds_x = [indsx(1:offsets(1,1)) indsx((offsets(1,1) + 2):2:(end-offsets(1,2))) indsx((end-offsets(1,2)+1):end)]; % psi-indices
            psiinds_y = [indsy(1:offsets(2,1)) indsy((offsets(2,1) + 2):2:(end-offsets(2,2))) indsy((end-offsets(2,2)+1):end)];
            phiinds_x = indsx((offsets(1,1) + 1):2:(end-offsets(1,2)));
            phiinds_y = indsy((offsets(2,1) + 1):2:(end-offsets(2,2)));
            
            resstart(1,res+1) = indsx(offsets(1,1)+1); resend(1,res+1)   = indsx(end-offsets(1,2));
            resstart(2,res+1) = indsy(offsets(2,1)+1); resend(2,res+1)   = indsy(end-offsets(2,2));
            
            sig_out(psiinds_x, phiinds_y,:) = sig_in( (endx-length(psiinds_x)+1):endx, 1:(endy-length(psiinds_y)), :);
            sig_out(phiinds_x, psiinds_y, :) = sig_in( 1:(endx-length(psiinds_x)), (endy-length(psiinds_y)+1):endy, :);
            sig_out(psiinds_x, psiinds_y, :) = sig_in( (endx-length(psiinds_x)+1):endx, (endy-length(psiinds_y)+1):endy, :);
            
            
            endx = endx - length(psiinds_x); endy = endy - length(psiinds_y);
            indsx = indsx((offsets(1,1)+1):2:(end-offsets(1,2))); 
            indsy = indsy((offsets(2,1)+1):2:(end-offsets(2,2)));
        end
        sig_out(indsx, indsy, :) = sig_in(1:endx, 1:endy, :);
    end
end
