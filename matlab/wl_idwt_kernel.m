function rekernel = wl_idwt_kernel(wave_name, data_size, dims, varargin)
    forward = 0;
    rekernel = wl_internal_kernel(wave_name, data_size, dims, forward, varargin{:});
end

