function x = wl_dwt2_impl_internal(x, fx, fy, prefilterx, prefiltery, offsets, m, bd_mode, data_layout)
    % Compute a 2D DWT using precomputed kernels. The kernels may be the default library kernels obtained by calling find_kernel, 
    % or may be used-defined.
    %
    % x:         Matrix whose DWT2 will be computed along the first dimensions. 
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
    
    indsx = 1:size(x,1); indsy = 1:size(x,2);

    % preconditioning   
    x(indsx, indsy, :) = tensor2_impl(x(indsx, indsy, :), @(x,bdm) prefilterx(x, 1), @(x,bdm) prefiltery(x, 1), bd_mode);
         
    for res = 0:(m - 1)
        x(indsx, indsy, :) = tensor2_impl(x(indsx, indsy, :), fx, fy, bd_mode);  
        indsx = indsx((offsets(1,1)+1):2:(end-offsets(1,2)));
        indsy = indsy((offsets(2,1)+1):2:(end-offsets(2,2)));
    end
    
    % postconditioning
    % x(indsx, indsy, :) = tensor2_impl(x(indsx, indsy, :), @(x,bdm) prefilterx(x, 0), @(x,bdm) prefiltery(x, 0), bd_mode);  
    
    x = reorganize_coeffs2_forward(x, m, offsets, data_layout);   
end     
        
function sig_out=reorganize_coeffs2_forward(sig_in, m, offsets, data_layout)
    sig_out = sig_in;
    if strcmpi(data_layout, 'resolution')
        indsx = 1:size(sig_in,1); indsy = 1:size(sig_in,2);
        endx = size(sig_in,1); endy = size(sig_in,2); 
        for res=1:m
            psiinds_x = [indsx(1:offsets(1,1)) indsx((offsets(1,1) + 2):2:(end-offsets(1,2))) indsx((end-offsets(1,2)+1):end)]; % psi-indices
            psiinds_y = [indsy(1:offsets(2,1)) indsy((offsets(2,1) + 2):2:(end-offsets(2,2))) indsy((end-offsets(2,2)+1):end)];
            phiinds_x = indsx((offsets(1,1) + 1):2:(end-offsets(1,2)));
            phiinds_y = indsy((offsets(2,1) + 1):2:(end-offsets(2,2)));
            
            sig_out( (endx-length(psiinds_x)+1):endx, 1:(endy-length(psiinds_y)), :) = sig_in(psiinds_x, phiinds_y,:);
            sig_out( 1:(endx-length(psiinds_x)), (endy-length(psiinds_y)+1):endy, :) = sig_in(phiinds_x, psiinds_y,:);
            sig_out( (endx-length(psiinds_x)+1):endx, (endy-length(psiinds_y)+1):endy, :) = sig_in(psiinds_x, psiinds_y, :);
            
            endx = endx - length(psiinds_x); endy = endy - length(psiinds_y);
            indsx = indsx((offsets(1,1)+1):2:(end-offsets(1,2)));
            indsy = indsy((offsets(2,1)+1):2:(end-offsets(2,2)));
        end
        sig_out(1:endx, 1:endy, :) = sig_in(indsx, indsy, :);
    end
end
