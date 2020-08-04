function x=wl_dwt_impl_from_kernel(x, rekernel, verbose)

    if nargin < 3
        verbose = 1;
    end

    if verbose
        x_size = size(x);

        if ~isequal(x_size, rekernel.data_size)
            disp('Error: wl_dwt_impl_from_kernel: Wrong input dimension.')    
        end

        if rekernel.forward ~=  1
            disp('Error: wl_dwt_impl_from_kernel: Not a forward kernel.')    
        end 

    end

    m = rekernel.m;
    bd_mode = rekernel.bd_mode;
    data_layout = rekernel.data_layout;

    if rekernel.dims == 1
        offsets = [rekernel.offset_Lx, rekernel.offset_Rx];
        if rekernel.transpose % if transpose, then f will we an idwt_kernel.
            x = wl_idwt1_impl_internal(x, rekernel.fx, rekernel.prefilterx, ...
                                       offsets, m, bd_mode, data_layout);
        else
            x =  wl_dwt1_impl_internal(x, rekernel.fx, rekernel.prefilterx, ...
                                       offsets, m, bd_mode, data_layout);
        end
    end

    if rekernel.dims == 2
        offsets = [rekernel.offset_Lx, rekernel.offset_Rx; ...
                   rekernel.offset_Ly, rekernel.offset_Ry];
        if rekernel.transpose % if transpose, then f will we an idwt_kernel, 
            x = wl_idwt2_impl_internal(x, rekernel.fx, rekernel.fy, ...
                   rekernel.prefilterx, rekernel.prefiltery, ...
                   offsets, m, bd_mode, data_layout);
        else
            x =  wl_dwt2_impl_internal(x, rekernel.fx, rekernel.fy, ...
                    rekernel.prefilterx, rekernel.prefiltery, ...
                    offsets, m, bd_mode, data_layout);
        end
    end

    if rekernel.dims == 3 % if not give error message
        offsets = [rekernel.offset_Lx, rekernel.offset_Rx; ...
                   rekernel.offset_Ly, rekernel.offset_Ry; ...
                   rekernel.offset_Lz, rekernel.offset_Rz];
        if rekernel.transpose % if transpose, then f will we an idwt_kernel, 
            x = wl_idwt3_impl_internal(x, rekernel.fx, rekernel.fy, ...
                   rekernel.fz, rekernel.prefilterx, ...
                   rekernel.prefiltery, rekernel.prefilterz, ...
                   offsets, m, bd_mode, data_layout);
        else
            x = wl_dwt3_impl_internal(x, rekernel.fx, rekernel.fy, ...
                   rekernel.fz, rekernel.prefilterx, ...
                   rekernel.prefiltery, rekernel.prefilterz, ...
                   offsets, m, bd_mode, data_layout);
        end
    end

end


