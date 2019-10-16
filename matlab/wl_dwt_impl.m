function x=wl_dwt_impl(x, wave_name, varargin)
    % Main function for computing the DWT of a given signal. Can be used for
    % all signals up to dimension 3.  The dimension of the data may be one
    % higher than the dimension of the transform, in which case the last
    % dimension is used for parallel computation.
    %
    % Note that this function computes all quantities needed from scratch in
    % order to compute the DWT for the wavelet in question.  This can be
    % time-consuming, and can be avoided by using the functions find_wav_props,
    % find_kernel,  the internal DWT functions dwt1_impl_internal,
    % dwt2_impl_internal, dwt3_impl_internal, as well as Matlabs persistence
    % functions. An example with minimum set of parameters is as follows:
    % 
    % [f, prefilter, offset_L, offset_R]=wl_find_kernel(wave_name, length_signal, 1)
    % save('wl_kernel.mat', 'f', 'prefilter', 'offset_L', 'offset_R');
    % ...
    % load('wl_kernel.mat');
    % x = wl_dwt1_impl_internal(x, f, prefilter, [offset_L offset_R]);
    %     
    % x:         Matrix whose DWT will be computed along the first dimension(s).      
    % wave_name: Name of the wavelet. Possible names are:
    %            'cdf97' - CDF 9/7 wavelet
    %            'spline53' - Spline 5/3 wavelet
    %            'splinex.x' - Spline wavelet with given number of vanishing 
    %                          moments for each filter
    %            'pwl0'  - Piecewise linear wavelet with 0 vanishing moments
    %            'pwl2'  - Piecewise linear wavelet with 2 vanishing moments
    %            'Haar'  - The Haar wavelet
    %            'dbX'   - Daubechies orthonormal wavelet with X vanishing
    %                      moments
    %            'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
    %                      with X vanishing moments
    %
    % This function also accepts a number of named, optional parameters. These are parsed by the function wl_setopts(). 
    % The documentation of this function also contains the full documentation of these optional parameters.

    args = {'wave_name', wave_name, 'dims', length(size(x))-1, varargin{:}};
    opts = wl_setopts(args{:});
                                                             
    [fx, prefilterx, offset_L, offset_R] = wl_find_kernel(wave_name, size(x,1), 1, args{:});
    offsets = [offset_L offset_R];
    if opts.dims == 1
        if opts.transpose % if transpose, then f will we an idwt_kernel.
            x = wl_idwt1_impl_internal(x, fx, prefilterx, offsets, args{:});
        else
            x =  wl_dwt1_impl_internal(x, fx, prefilterx, offsets, args{:});
        end
    else
        [fy, prefiltery, offset_L, offset_R] = wl_find_kernel(wave_name, size(x,2), 1, args{:});
        offsets = [offsets; offset_L offset_R];
        if opts.dims == 2
            if opts.transpose % if transpose, then f will we an idwt_kernel, 
                x = wl_idwt2_impl_internal(x, fx, fy, prefilterx, prefiltery, offsets, args{:});
            else
                x =  wl_dwt2_impl_internal(x, fx, fy, prefilterx, prefiltery, offsets, args{:});
            end
        else
            [fz, prefilterz, offset_L, offset_R] = wl_find_kernel(wave_name, size(x,3), 1, args{:});
            offsets = [offsets; offset_L offset_R];
            if opts.dims == 3 % if not give error message
                if opts.transpose % if transpose, then f will we an idwt_kernel, 
                    x = wl_idwt3_impl_internal(x, fx, fy, fz, prefilterx, prefiltery, prefilterz, offsets, args{:});
                else
                    x =  wl_dwt3_impl_internal(x, fx, fy, fz, prefilterx, prefiltery, prefilterz, offsets, args{:});
                end
            end
        end
    end         
end
