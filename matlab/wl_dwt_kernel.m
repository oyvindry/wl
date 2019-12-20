function rekernel = wl_dwt_kernel(wave_name, data_size, dims, varargin)
    forward = 1;
    rekernel = wl_internal_kernel(wave_name, data_size, dims, forward, varargin{:});
end

