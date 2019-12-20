function rekernel = wl_internal_kernel(wave_name, data_size, dims, forward, varargin)
    args = {'wave_name', wave_name, 'dims', dims, varargin{:}};
    opts = wl_setopts(args{:});

    if iscolumn(data_size)
        data_size = data_size';
    end

    rekernel = opts;
    rekernel.data_size = data_size;
    rekernel.forward = forward;
    rekernel.args = args;

    if opts.dims >= 1    
        [fx, prefilterx, offset_L, offset_R] = wl_find_kernel(wave_name, data_size(1), forward, args{:});

        rekernel.fx = fx;
        rekernel.prefilterx = prefilterx;
        rekernel.offset_Lx = offset_L;
        rekernel.offset_Rx = offset_R;
    end

    if opts.dims >= 2
        [fy, prefiltery, offset_L, offset_R] = wl_find_kernel(wave_name, data_size(2), forward, args{:});

        rekernel.fy = fy;
        rekernel.prefiltery = prefiltery;
        rekernel.offset_Ly = offset_L;
        rekernel.offset_Ry = offset_R;
    end
    
    if opts.dims >= 3
        [fz, prefilterz, offset_L, offset_R] = wl_find_kernel(wave_name, data_size(3), forward, args{:});

        rekernel.fz = fz;
        rekernel.prefilterz = prefilterz;
        rekernel.offset_Lz = offset_L;
        rekernel.offset_Rz = offset_R;
    end
end
