function opts=wl_setopts(varargin)
    % Set defaul parameters
    
    % m:              Number of resolutions. Default: 1.
    % bd_mode:        Boundary extension mode. Possible modes are. 
    %                 'per'    - Periodic extension (default for wavelets with non-symmetric filters)
    %                 'symm'   - Symmetric extension (default for wavelets with symmetric filters)
    %                 'none'   - Take no extra action at the boundaries (i.e., the boundaries are zero-padded)
    %                 'bd'     - Preserve vanishing moments at the boundaries
    %                 If this parameter is not provided, the implementation looks for a given wavelet name, and uses the default bd_mode for that wavelet.  
    %                 If no wavelet is given, the implementation looks for given filters, and uses 'symm' or 'per' depending on whether the filters are symmetric or not.
    %                 If no filters are given, the 'symm' mode is chosen.  
    % prefilter_mode: Possible modes are:
    %                 'none' (default)
    %                 'filter'
    %                 'bd_pre' - Preconditioning. Applies only in combination with preservation of vanishing moments at the boundaries.
    % dims:           The number of dimensions to apply the DWT to. Always applied
    %                 to the first dimensions. Default: max(numdims(x)-1,1). This means 
    %                 that sound with many channels, and images with many colour 
    %                 components default to a one- and two-dimensional DWT, 
    %                 respectively
    % dual:           Whether to apply the dual wavelet rather than the wavelet 
    %                 itself. Default: 0
    % transpose:      Whether the transpose is to be taken. Default: 0
    % data_layout:    How data should be assembled. Possible modes are:
    %                 'resolution': Lowest resolution first (default)
    %                 'time': Sort according to time
    % impl_strategy:  How the kernel should be implemented. Possible modes are:
    %                 'lifting'. Apply a lifting-based procedure. An error is returned if the filters can't be factored in terms of elementary lifting steps.
    %                 'filter'. Apply a filter-based procedure. This mode works for any wavelet.
    %                 'any' (default). With this setting the library will find the optimal implementation on its own. 
    %                 As a rule of thumb, a lifting-based implementation is applied for wavelets with least dissimilar filters. 
    %                 Otherwise a filter-based implementation is applied.
    % symbolic        Whether symbolic computation should be attempted when computing the boundary coefficients. Defaults to 1 for spline wavelets. 
    %                 Setting this to 1 will increase the precomputation time considerably. 
    % staggered       Whether staggered supports should be made
    % polbasis        Which polynomial basis to use. Bernstein, Gramm.
    
    % Parse input arguments. Set default parameters
    opts.dims = 0;          
    opts.m = 1;
    opts.prefilter_mode = 'none';
    opts.dual = 0;
    opts.transpose = 0;
    opts.data_layout = 'resolution';
    opts.impl_strategy = 'any';
    opts.wave_name = 'unknown';
    opts.symbolic = 0;
    opts.staggered = 1;
    opts.polbasis = 'gramm';
    
    args = varargin;
    if ~iscell(args)
        args = {args};
    end
    
    nbr_args = numel(args);
    for i = 1:2:nbr_args
        if (ischar(args{i}))
            if strcmpi(args{i},'m')
                if i+1 <= nbr_args
                    if isnumeric(args{i+1})
                        opts.m = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'bd_mode')
                if i+1 <= nbr_args
                    if ischar(args{i+1})
                        opts.bd_mode = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'prefilter_mode')
                if i+1 <= nbr_args
                    if ischar(args{i+1})
                        opts.prefilter_mode = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'dual')
                if i+1 <= nbr_args
                    if isnumeric(args{i+1})
                        opts.dual = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'transpose')
                if i+1 <= nbr_args
                    if isnumeric(args{i+1})
                        opts.transpose = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'dims')
                if i+1 <= nbr_args
                    if isnumeric(args{i+1})
                        if args{i+1} ~= 0
                            opts.dims = args{i+1};    
                        end
                    end
                end
            end
            if strcmpi(args{i},'data_layout')
                if i+1 <= nbr_args
                    if ischar(args{i+1})
                        opts.data_layout = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'impl_strategy')
                if i+1 <= nbr_args
                    if ischar(args{i+1})
                        opts.impl_strategy = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'wave_name')
                if i+1 <= nbr_args
                    if ischar(args{i+1})
                        opts.wave_name = args{i+1};
                        if strcmpi(opts.wave_name(1:2), 'db') | strcmpi(opts.wave_name, 'Haar') | strcmpi(opts.wave_name(1:3), 'sym')
                            opts.bd_mode = 'per';
                        else 
                            opts.bd_mode = 'symm';
                        end
                    end
                end
            end
            if strcmpi(args{i},'symbolic')
                if i+1 <= nbr_args
                    if isnumeric(args{i+1})
                        opts.symbolic = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'staggered')
                if i+1 <= nbr_args
                    if isnumeric(args{i+1})
                        opts.staggered = args{i+1};
                    end
                end
            end
            if strcmpi(args{i},'polbasis')
                if i+1 <= nbr_args
                    if ischar(args{i+1})
                        opts.polbasis = args{i+1};
                    end
                end
            end
        end
    end
    % opts % print arguments
end
