function [f, prefilter, offset_L, offset_R]=find_kernel_from_filters(g0, N, h0, Ntilde, length_signal, forward, varargin)
    % Function which returns the default kernel of the library for use with the wavelet with properties encapsulated by the given parameters.
    % The kernel can be passed to the the internal functions (i)dwt1_impl_internal, (i)dwt2_impl_internal, (i)dwt3_impl_internal. 
    % User-defined kernels can also be passed to these internal functions: They simply have to take the x and bd_mode parameters, and return the 
    % transformed vectors.
    %
    % It is assumed that the wavelet is delay-normalized, with alpha=-1 (even length filters) or alpha=1 (odd length filters). 
    % These values of alpha make the transforms compatible with standard such transforms as they are defined in the literature. These requirements also define the high-pass filters uniquely from the given low-pass filters.
    %
    % The assumptions on the center of the supports removes the need to indicate the position of the filter coefficient with index 0. If the support centers for (g0,h0) are (0,0) or (1/2,-1/2), 
    % the boundary wavelet construction is hugely simplified, as in the paper. 
    %
    % g0:             Low pass synthesis filter. It is assumed that the support of the filter coefficients is symmetric about 0 or 1/2, depending on whether there are an even or odd number of filter coefficients.
    % N:              Number of vanishing moments
    % h0:             Low pass analysis filter. It is assumed that the support of the filter coefficients is symmetric about 0 or -1/2, depending on whether there are an even or odd number of filter coefficients.
    % Ntilde:         Number of vanishing moments in the dual transform
    % length_signal:  Length of the input signal.
    % forward:        Whether the forward transform should be used
    %
    % This function also accepts a number of named, optional parameters. These are parsed by the function wl_setopts(). 
    % The documentation of this function also contains the full documentation of these optional parameters.
    
    if mod(length(h0),2)==1 % TODO: This test should be on symmetry rather than filter lengths
        bd_mode = 'symm';
    else
        bd_mode = 'per';
    end
    opts = wl_setopts('wave_name', 'unknown', 'bd_mode', bd_mode, varargin{:});
    
    if mod(length(g0),2)==0
        %g1=h0.*(-1).^(0:(length(g0)-1));
        %h1=g0.*(-1).^(1:(length(g0)));
        g1=-h0.*(-1).^((-length(g0)/2):(length(g0)/2-1));
        h1=-g0.*(-1).^((-(length(g0)/2-1)):(length(g0)/2)); % alpha=-1
    else
        h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2)); % alpha=1
        g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
    end
    
    if strcmpi(opts.impl_strategy,'any')
        if abs(length(g0)-length(g1))<=2 % Least dissimilar
            opts.impl_strategy = 'lifting';
        else
            opts.impl_strategy = 'filter';
        end
    end
    
    prefilter = @(x, forward) x;
    if opts.transpose
        forward = ~forward;
        opts.dual = ~opts.dual;
    end
    if opts.dual 
        [dual_wav_props, wav_props, WLtilde, WL, WRtilde, WR]=find_wav_props_from_filters(g0, g1, N, h0, h1, Ntilde, opts.m, opts.bd_mode, length_signal, opts.impl_strategy, opts.prefilter_mode);
    else
        [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde]=find_wav_props_from_filters(g0, g1, N, h0, h1, Ntilde, opts.m, opts.bd_mode, length_signal, opts.impl_strategy, opts.prefilter_mode);
    end
    offset_L = wav_props.offset_L;
    offset_R = wav_props.offset_R;
    
    if strcmpi(opts.bd_mode, 'bd') % tail handling matrices have not been computed yet
        if forward
            if strcmpi(opts.impl_strategy, 'lifting')
                fnew = @(x, bdm) idwt_kernel_lifting(x, bdm, dual_wav_props);
            elseif strcmpi(opts.impl_strategy, 'filter')
                fnew = @(x, bdm) idwt_kernel_filters(x, bdm, dual_wav_props);
            end
            [dual_wav_props.A_L,dual_wav_props.A_R]=find_AL_AR(WLtilde, WRtilde, dual_wav_props, fnew); % fnew is Gtilde
            f = find_kernel_general(opts.prefilter_mode, opts.impl_strategy, @(x, bdm) dwt_kernel_lifting(x, bdm, dual_wav_props), @(x, bdm) dwt_kernel_filters(x, bdm, dual_wav_props));
        else
            if strcmpi(opts.impl_strategy, 'lifting')
                fnew = @(x, bdm) idwt_kernel_lifting(x, bdm, wav_props);
            elseif strcmpi(opts.impl_strategy, 'filter')
                fnew = @(x, bdm) idwt_kernel_filters(x, bdm, wav_props);
            end
            [wav_props.A_L,wav_props.A_R]=find_AL_AR(WL, WR, wav_props, fnew); % fnew is G
            f = find_kernel_general(opts.prefilter_mode, opts.impl_strategy, @(x, bdm) idwt_kernel_lifting(x, bdm, wav_props), @(x, bdm) idwt_kernel_filters(x, bdm, wav_props));
        end
        if strcmpi(opts.prefilter_mode, 'bd_pre') 
            prefilter = @(x, fwd) precond_impl(x, fwd, wav_props);
        elseif strcmpi(opts.prefilter_mode, 'filter')
            [wav_props.filtermatr, wav_props.d] = compute_filter_matrix(wav_props, N, max(N,Ntilde), length_signal);
            prefilter = @(x, fwd) prefilter_impl(x, fwd, wav_props);
        elseif ~strcmpi(opts.prefilter_mode, 'none')
            throw(MException('WL:illegal_parameters', sprintf('Illegal prefiltering mode: %s', opts.prefilter_mode)));
        end
    else
        if forward
            f = find_kernel_general(opts.prefilter_mode, opts.impl_strategy, @(x, bdm) dwt_kernel_lifting(x, bdm, dual_wav_props), @(x, bdm) dwt_kernel_filters(x, bdm, dual_wav_props));
        else
            f = find_kernel_general(opts.prefilter_mode, opts.impl_strategy, @(x, bdm) idwt_kernel_lifting(x, bdm, wav_props), @(x, bdm) idwt_kernel_filters(x, bdm, wav_props));
        end
    end
end

function [A_L,A_R]=find_AL_AR(WL, WR, wp, f)
    WL = double(WL); WR= double(WR);
    M = max(size(WL)+size(WR));
    if mod(M-wp.length_signal,2) == 1
        M = M + 1;
    end
    
    x1 = eye(M);
    x1 = wl_idwt1_impl_internal(x1, f, @(x, bdm) x, [wp.offset_L wp.offset_R], 'bd_mode', 'none', 'data_layout', 'time');
    [w1, w2] = size(WL); A_L = WL - x1(1:w1,1:w2);
    [w1, w2] = size(WR); A_R = WR - x1((M-w1+1):M,(M-w2+1):M);
    A_L = double(A_L); A_R = double(A_R);
end


function f = find_kernel_general(prefilter_mode, impl_strategy, kernel_lifting, kernel_filters)
    f = @(x, bd_mode) x;
    if strcmpi(impl_strategy, 'lifting')
        f = kernel_lifting;
    elseif strcmpi(impl_strategy, 'filter')
        f = kernel_filters;
    end
end

function x=filter_impl(t, x, bd_mode, swap)
    if mod(length(t),2)==0
        if swap
            l = length(t)/2; r = l-1;
        else
            r = length(t)/2; l = r-1;
        end
    else
        r = (length(t) - 1)/2; l = r;
    end
    
    N = size(x, 1);
    szx = size(x);
    n = prod(szx(2:end));
    
    if strcmpi(bd_mode, 'symm')
        assert(l==r);
        y = [x((l+1):(-1):2, :) ; x(:,:); x((N-1):(-1):(N - l), :)];
    elseif strcmpi(bd_mode, 'per')
        y = [x((N - r + 1):N, :); x(:,:); x(1:l, :)];
    elseif strcmpi(bd_mode, 'none') || strcmpi(bd_mode, 'bd')
        y = [ zeros(r, n); x(:, :); zeros(l, n)];
    end
    for k=1:size(y,2)
        z = conv(t, y(:,k));
        x(:,k) = z(length(t):(length(z)-length(t)+1));
    end
end

function x=dwt_kernel_filters(x, bd_mode, wp)
    if strcmpi(bd_mode, 'bd')
        y1 = wp.A_L'*x(1:size(wp.A_L,1), :);
        y2 = wp.A_R'*x((end-size(wp.A_R,1)+1):end, :);
    end
    x0 = filter_impl(flip(wp.g0), x, bd_mode, ~wp.swap);
    x1 = filter_impl(flip(wp.g1), x, bd_mode, wp.swap);
    x(1:2:end,:) = x0(1:2:end,:);
    x(2:2:end,:) = x1(2:2:end,:);
    if strcmpi(bd_mode, 'bd')
        x(1:size(wp.A_L,2), :) ...
            = x(1:size(wp.A_L,2), :) + y1;
        x((end-size(wp.A_R,2)+1):end, :) ...
            = x((end-size(wp.A_R,2)+1):end, :) + y2;
    end
end
% End dwt_kernel_filters
    
function x=idwt_kernel_filters(x, bd_mode, wp)
    if strcmpi(bd_mode, 'bd')
        y1 = wp.A_L*x(1:size(wp.A_L,2), :);
        y2 = wp.A_R*x((end-size(wp.A_R,2)+1):end, :);
    end
    x0 = x; x0(2:2:end,:) = 0;
    x1 = x; x1(1:2:end,:) = 0;
    x0 = filter_impl(wp.g0, x0, bd_mode, wp.swap);
    x1 = filter_impl(wp.g1, x1, bd_mode, ~wp.swap);
    x = x0 + x1;
    if strcmpi(bd_mode, 'bd')
        x(1:size(wp.A_L,1), :) ...
            = x(1:size(wp.A_L,1), :) + y1;
        x((end-size(wp.A_R,1)+1):end, :) ....
            = x((end-size(wp.A_R,1)+1):end, :) + y2;
    end
end
% End idwt_kernel_filters

function x=dwt_kernel_lifting(x, bd_mode, wp)
    if wp.lambdas(:,1) == wp.lambdas(:,2)
        f_even = @(lambda1, lambda2, x, bd_mode) lifting_even_symm(lambda1, x, bd_mode);
        f_odd =  @(lambda1, lambda2, x, bd_mode) lifting_odd_symm(lambda1, x, bd_mode);
    else
        f_even = @lifting_even;
        f_odd = @lifting_odd;
    end
    if strcmpi(bd_mode, 'bd')
        y1 = wp.A_L'*x(1:size(wp.A_L,1), :);
        y2 = wp.A_R'*x((end-size(wp.A_R,1)+1):end, :);
    end
    x(1:2:end, :) = x(1:2:end, :)/wp.alpha;
    x(2:2:end, :) = x(2:2:end, :)/wp.beta;
    iseven = ~wp.last_even;
    for stepnr = (size(wp.lambdas, 1)):(-1):1
        if iseven
            x = f_even(wp.lambdas(stepnr,2), wp.lambdas(stepnr,1), x, bd_mode);
        else
            x = f_odd(wp.lambdas(stepnr,2), wp.lambdas(stepnr,1), x, bd_mode);
        end
        iseven = ~iseven;
    end
    if strcmpi(bd_mode, 'bd')
        x(1:size(wp.A_L,2), :) = x(1:size(wp.A_L,2), :) + y1;
        x((end-size(wp.A_R,2)+1):end, :) = x((end-size(wp.A_R,2)+1):end, :) + y2;
    end
end
% End dwt_kernel_lifting

function x=idwt_kernel_lifting(x, bd_mode, wp)
    if wp.lambdas(:,1) == wp.lambdas(:,2) 
        f_even = @(lambda1, lambda2, x, bd_mode) lifting_even_symm(lambda1, x, bd_mode);
        f_odd =  @(lambda1, lambda2, x, bd_mode) lifting_odd_symm(lambda1, x, bd_mode);
    else
        f_even = @lifting_even;
        f_odd = @lifting_odd;
    end
    if strcmpi(bd_mode, 'bd')
       y1 = wp.A_L*x(1:size(wp.A_L,2), :);
       y2 = wp.A_R*x((end-size(wp.A_R,2)+1):end, :);
    end
    iseven = ( mod(size(wp.lambdas, 1), 2) == wp.last_even );
    for stepnr = 1:(size(wp.lambdas, 1))
        if iseven
            x = f_even(wp.lambdas(stepnr, 1), wp.lambdas(stepnr, 2), x, bd_mode);
        else
            x = f_odd(wp.lambdas(stepnr, 1), wp.lambdas(stepnr, 2), x, bd_mode);
        end
        iseven = ~iseven;
    end
    x(1:2:end, :)=x(1:2:end, :)/wp.alpha;
    x(2:2:end, :)=x(2:2:end, :)/wp.beta;

    if strcmpi(bd_mode, 'bd')
        x(1:size(wp.A_L,1), :) = ...
            x(1:size(wp.A_L,1), :) + y1;
        x((end-size(wp.A_R,1)+1):end, :) = ...
            x((end-size(wp.A_R,1)+1):end, :) + y2;
    end
end
% End idwt_kernel_lifting

function x=precond_impl(x, forward, wp)
    n = size(wp.A_L_pre_inv,1);
    if forward == 1
        x(1:n,:)           = wp.A_L_pre_inv\x(1:n,:);
        x((end-n+1):end,:) = wp.A_R_pre_inv\x((end-n+1):end,:);
    else
        x(1:n,:)           = wp.A_L_pre_inv*x(1:n,:);
        x((end-n+1):end,:) = wp.A_R_pre_inv*x((end-n+1):end,:);
    end
end

function x=prefilter_impl(x, forward, wp) 
   if forward == 1
        x=solve_A(wp.filtermatr, x, wp.d);
    else
        x=multiply_A(wp.filtermatr, x, wp.d);
    end
end

% Assumes that A is obtained from luimpl_banded(A, d)
function x=solve_A(A, b, d)
    z=rforwardsolve(A, b, d);
    x=rbacksolve(A, z, d);
end

function b=multiply_A(A, x, d)
    n = size(x,1);
    z = x;
    for k=1:n
        currentrow = mod(k - 1, 2*d + 1) + 1;
        uk=min(n,k+d);
        z(k,:) = A(currentrow, k:uk)*x(k:uk,:);
    end
    b = z;
    for k=2:n
        currentrow = mod(k - 1, 2*d + 1) + 1;
        lk = max(1,k-d);
        b(k,:) = A(currentrow, lk:(k-1))*z(lk:(k-1),:) + z(k,:);
    end
end

function x=rforwardsolve(A, b, d)
    n=size(b,1); x=b;
    x(1,:)=b(1,:);
    for k=2:n
        currentrow = mod(k - 1, 2*d + 1) + 1;
        lk=max(1,k-d);
        x(k,:)=b(k,:)-A(currentrow,lk:(k-1))*x(lk:(k-1),:);
    end
end
   
function x=rbacksolve(A, b, d)
    n=size(b, 1); x=b;
    x(n,:)=b(n,:)/A(mod(n - 1, 2*d + 1) + 1,n);
    for k=n-1:-1:1
        currentrow = mod(k - 1, 2*d + 1) + 1;
        uk=min(n,k+d);
        x(k,:)=(b(k,:)-A(currentrow,(k+1):uk)*x((k+1):uk,:))/A(currentrow,k);
    end
end




function x=lifting_even_symm(lambda1, x, bd_mode)
    N = size(x, 1);
    if strcmpi(bd_mode, 'per')
        assert(mod(N,2) == 0)
    end
    if strcmpi(bd_mode, 'symm')
        x(1, :) = x(1, :) + 2*lambda1*x(2, :); % Symmetric extension
    elseif strcmpi(bd_mode, 'per')
        x(1, :) = lambda1*(x(2, :) + x(N, :)) + x(1, :);
    elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
        x(1, :) = lambda1*x(2, :) + x(1, :);
    end
    x(3:2:(N-1), :) = x(3:2:(N-1),:) + lambda1*(x(2:2:(N-2),:) + x(4:2:N,:));
    if mod(N,2) == 1 % last must also be included
        if strcmpi(bd_mode, 'symm')
            x(N, :) = x(N, :) + 2*lambda1*x(N-1, :); % Symmetric extension
        elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
            x(N, :) = x(N, :) + lambda1*x(N-1, :);
        end
    end
end
% End lifting_even_symm
    
function x=lifting_odd_symm(lambda1, x, bd_mode)
    N = size(x, 1);
    x(2:2:(N-1), :) = x(2:2:(N-1),:) + lambda1*(x(1:2:(N-2),:) + x(3:2:N,:));
    if mod(N,2)==0 % last must also be included
        if strcmpi(bd_mode, 'symm')
            x(N, :) = x(N, :) + 2*lambda1*x(N-1, :); % Symmetric extension
        elseif strcmpi(bd_mode, 'per')
            x(N, :) = lambda1*(x(1, :) + x(N-1, :)) + x(N, :);
        elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
            x(N, :) = lambda1*x(N-1, :) + x(N, :);
        end
    end
end
% End lifting_odd_symm

function x=lifting_even(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    if strcmpi(bd_mode, 'per')
        x(1, :) = lambda1*x(2, :) + x(1, :) + lambda2*x(N, :);
    elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
        x(1, :) = lambda1*x(2, :) + x(1, :);
    end
    x(3:2:(N-1), :) = ...
        lambda1*x(4:2:N, :) + x(3:2:(N-1), :) + lambda2*x(2:2:(N-2), :);
    if strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
        if mod(N,2) == 1 % last must also be included
            x(N, :) = x(N, :) + lambda2*x(N-1, :);
        end
    end
end
% End lifting_even

function x=lifting_odd(lambda1, lambda2, x, bd_mode)
    N = size(x, 1);
    x(2:2:(N-1), :) = ...
        lambda1*x(3:2:N, :) + x(2:2:(N-1), :) + lambda2*x(1:2:(N-2), :);
    if mod(N,2) == 0 % last must also be included
        if strcmpi(bd_mode, 'per')
            x(N, :) = lambda1*x(1, :) + x(N, :) + lambda2*x(N-1, :);
        elseif strcmpi(bd_mode, 'bd') || strcmpi(bd_mode, 'none')
            x(N, :) = x(N, :) + lambda2*x(N-1, :);
        end
    end
end
% End lifting_odd






function [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde]=find_wav_props_from_filters(g0, g1, N, h0, h1, Ntilde, m, bd_mode, length_signal, impl_strategy, prefilter_mode)
    % Computes the properties of a wavelet with a given set of filters. What properties 
    % are computed depend on the bd_mode parameter, m, and length_signal.
    %
    % g0:        Low-pass filter coefficients
    % g1:        High-pass filter coefficients
    % m:         Number of resolutions. Default: 1
    % bd_mode:   Boundary extension mode. Possible modes are. 
    %            'per'    - Periodic extension
    %            'symm'   - Symmetric extension (default)
    %            'none'   - Take no extra action at the boundaries
    %            'bd'     - Boundary wavelets
    % length_signal: Length of the input signal. Default: 0.
    % impl_strategy: lifting, filter, any (default)
     
    wav_props.m = m; dual_wav_props.m = m;
    wav_props.length_signal = length_signal; dual_wav_props.length_signal = length_signal;
    wav_props.offset_L = 0; dual_wav_props.offset_L = 0;
    wav_props.offset_R = 0; dual_wav_props.offset_R = 0;
    wav_props.g0 = g0; wav_props.g1 = g1;
    dual_wav_props.g0 = flip(h0); dual_wav_props.g1 =flip(h1);
    wav_props.N = N; dual_wav_props.N = Ntilde;
    wav_props.swap = 0; dual_wav_props.swap = 0;
    
    % Compute values of scaling function at the integers
    [wav_props.L,wav_props.R]=findsupports(wav_props.g0);
    inds = (wav_props.L+1):(wav_props.R-1);
    matr=Gsegment(wav_props.g0, wav_props.L:wav_props.R, 2*inds, inds, 0);
    eigval = sum(wav_props.g0)/2;
    wav_props.phivals = null(matr-eigval*eye(size(matr)));
    % assert(size(wav_props.phivals,2)==1);
    wav_props.phivals = wav_props.phivals/sum(wav_props.phivals);
    % wav_props.phivals
    
    % Compute values of dual scaling function at the integers
    [dual_wav_props.L,dual_wav_props.R]=findsupports(dual_wav_props.g0);
    inds = (dual_wav_props.L+1):(dual_wav_props.R-1);
    matr=Gsegment(dual_wav_props.g0, dual_wav_props.L:dual_wav_props.R, 2*inds, inds, 0);
    eigval = sum(dual_wav_props.g0)/2;
    dual_wav_props.phivals = null(matr-eigval*eye(size(matr)));
    % assert(size(dual_wav_props.phivals,2)==1); %TODO: This assert kicks in for pwl0
    dual_wav_props.phivals = dual_wav_props.phivals/sum(dual_wav_props.phivals);
    % dual_wav_props.phivals
     
    % TODO: check perfect reconstruction
    
    WL = 0; WLtilde = 0; WR = 0; WRtilde = 0;
    if strcmpi(bd_mode, 'bd')
        [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde] = wav_props_general(N, Ntilde, wav_props, dual_wav_props, prefilter_mode); 
    end 
    
    if strcmpi(impl_strategy,'lifting')
        [wav_props.lambdas, wav_props.alpha, wav_props.beta, wav_props.last_even]=lifting_fact(h0, h1);
        dual_wav_props.lambdas = -fliplr(wav_props.lambdas);
        dual_wav_props.alpha = 1/wav_props.alpha;
        dual_wav_props.beta = 1/wav_props.beta;
        dual_wav_props.last_even = ~wav_props.last_even;
        if strcmpi(bd_mode, 'bd')
            if mod(wav_props.offset_L, 2) == 1
                wav_props.last_even = ~wav_props.last_even;
                betaval = wav_props.alpha; wav_props.alpha = wav_props.beta; wav_props.beta = betaval;
                dual_wav_props.last_even = ~dual_wav_props.last_even;
                betaval = dual_wav_props.alpha; dual_wav_props.alpha = dual_wav_props.beta; dual_wav_props.beta = betaval;
            end
        end 
    end
    
    if strcmpi(bd_mode, 'bd')
        % swap filters if offset is odd
        if mod(wav_props.offset_L,2) == 1
            wav_props.swap = 1; dual_wav_props.swap = 1;
            g0temp = wav_props.g0; wav_props.g0 = wav_props.g1; wav_props.g1 = g0temp;
            g0temp = dual_wav_props.g0; dual_wav_props.g0 = dual_wav_props.g1; dual_wav_props.g1 = g0temp;
        end
    end
end

function [lambdas, alpha, beta, last_even] = lifting_fact(h0, h1)
    % Compute a lifting factorization so that
    % [alpha 0; beta 0] = \Lambda_n \cdots \Lambda_1 H
    % lambdas: [\Lambda_1; \Lambda_2; \cdots \Lambda_n]
    % last_even: If \Lambda_n is even
    
    lambdas = [];
    if mod(length(h0),2)==0 % L+R=1
        if length(h0)~=length(h1)
            if mod(length(h0)/2,2)==0, 
                throw(MException('WL:lifting_impossible', 'Cant make lifting factorization'));
            end
        end
        
        if mod(length(h0)/2, 2)==0
            h00=h0(1:2:end); h01=h0(2:2:end); 
        else 
            h00=h0(2:2:end); h01=h0(1:2:end);
        end
        
        if mod(length(h1)/2, 2)==0
            h10=h1(1:2:end); h11=h1(2:2:end);
        else 
            h10=h1(2:2:end); h11=h1(1:2:end);
        end
        
        if length(h0) == length(h1)
            if mod(length(h0)/2, 2)==0
                lambda1=-h00(1)/h10(1); lambda2 = 0;
                h00=h00+lambda1*h10; h00 = h00(2:end);
                h01=h01+lambda1*h11; h01 = h01(2:end);
            else
                lambda1 = 0; lambda2=-h10(end)/h00(end); 
                h10=h10+lambda2*h00; h10 = h10(1:(end-1));
                h11=h11+lambda2*h01; 
                if length(h11)>1 
                    h11 = h11(1:(end-1));
                end
            end
            lambdas=[lambdas; [lambda1 lambda2]];
        end
    else % L+R=0
        if mod((length(h0)-1)/2, 2)==0
            h00=h0(1:2:end); h01=h0(2:2:end);
        else
            h00=h0(2:2:end); h01=h0(1:2:end);
        end
        
        if mod((length(h1)-1)/2, 2)==0
            h10=h1(2:2:end); h11=h1(1:2:end);
        else
            h10=h1(1:2:end); h11=h1(2:2:end);
        end
    end
    
    if abs(length(h00)-length(h10))~=1
        throw(MException('WL:lifting_impossible', 'Wavelet not least dissimilar'));
    end
    if abs(length(h01)-length(h11))>1
        throw(MException('WL:lifting_impossible', 'Wavelet not least dissimilar'));
    end
    
    while length(h10)>0 & length(h01)>0
        if length(h00)>length(h10) % Reduce the degree in the first row.
            lambda1=-h00(1)/h10(1); lambda2=-h00(end)/h10(end);
            h00=h00+conv(h10,[lambda1 lambda2]); h00 = h00(2:(end-1));
            h01=h01+conv(h11,[lambda1 lambda2]); h01 = h01(2:(end-1));
        else % reduce the degree in the second row. 
            lambda1=-h10(1)/h00(1); lambda2=-h10(end)/h00(end);
            h10=h10+conv(h00,[lambda1 lambda2]); h10 = h10(2:(end-1));
            h11=h11+conv(h01,[lambda1 lambda2]); 
            if length(h11)==2
                h11=sum(h11);
            else
                h11 = h11(2:(end-1));
            end
        end
        lambdas=[lambdas; [lambda1 lambda2]];
    end
  
    % Add the final lifting, and alpha,beta
    alpha=h00;
    beta=h11;
    if length(h10)==0 % zero in the lower left
        lastlift = -h01/beta;
        last_even = 1;
    else % zero in the upper right
        lastlift = -h10/alpha;
        last_even = 0;
    end
    if mod(length(h0),2)==0
        if mod(length(h0)/2,2)==0
            lambdas = [lambdas; [0 lastlift]];
        else
            lambdas = [lambdas; [lastlift 0]];
        end
    else
        lambdas = [lambdas; lastlift];
    end
end

function [L,R]=findsupports(g0)
     % Find L, R, from g0
    if mod(length(g0),2) == 0 
        k = length(g0)/2;
        L = -k+1; R = k;
    else
        k = (length(g0)-1)/2;
        L = -k; R = k;
    end
end

function [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde]=wav_props_general(N, Ntilde, wav_props, dual_wav_props, prefilter_mode)
    Nprime = max(N, Ntilde);
    
    % Define K and Ktilde so that they satisfy the requirements of Definition 3.1.
    wav_props.K_L = max(-wav_props.L,N); 
    dual_wav_props.K_L = max(-dual_wav_props.L,Ntilde);
    if wav_props.K_L-N < dual_wav_props.K_L - Ntilde
        wav_props.K_L = dual_wav_props.K_L - Ntilde + N;
    end
    if wav_props.K_L-N > dual_wav_props.K_L-Ntilde
        dual_wav_props.K_L = wav_props.K_L-N+Ntilde;
    end
    wav_props.K_R      = wav_props.K_L;
    dual_wav_props.K_R = dual_wav_props.K_L;
    
    % Make sure that K and Ktilde satisfy the requirement of Theorem 6.2, and are as equal as possible.
    s = wav_props.length_signal + wav_props.R + wav_props.L - 2*N - 1 + wav_props.K_L + wav_props.K_R;
    if mod(s,2^wav_props.m) > 0
        toadd = 2^wav_props.m - mod(s,2^wav_props.m);
        wav_props.K_L = wav_props.K_L + floor(toadd/2); dual_wav_props.K_L = dual_wav_props.K_L + floor(toadd/2);
        wav_props.K_R = wav_props.K_R + ceil(toadd/2);  dual_wav_props.K_R = dual_wav_props.K_R + ceil(toadd/2);
    end
    
    % Set offsets
    wav_props.offset_L = wav_props.K_L - N; wav_props.offset_R = wav_props.K_R - N; dual_wav_props.offset_L = wav_props.offset_L; dual_wav_props.offset_R = wav_props.offset_R;
    
    % The left edge
    [WL, WLtilde, wav_props.A_L_pre_inv, dual_wav_props.A_L_pre_inv, N0L, wav_props.CL, dual_wav_props.CL] = bw_compute_left(wav_props.g0, wav_props.g1, N, wav_props.K_L, dual_wav_props.g0, dual_wav_props.g1, Ntilde, dual_wav_props.K_L);
    
    % The right edge
    [WR, WRtilde, wav_props.A_R_pre_inv, dual_wav_props.A_R_pre_inv, N0R, wav_props.CR, dual_wav_props.CR] = bw_compute_left(flip(wav_props.g0), flip(wav_props.g1), N, wav_props.K_R, flip(dual_wav_props.g0), flip(dual_wav_props.g1), Ntilde, dual_wav_props.K_R);

    dimphi1 = 2^(1-wav_props.m)*wav_props.length_signal + (1-2^(1-wav_props.m))*(2*N-wav_props.L-wav_props.R-wav_props.K_L-wav_props.K_R+1); % Equation (6.4)
    wav_props.s_L = wav_props.K_L - N + 2*max(Nprime,N0L); dual_wav_props.s_L = wav_props.s_L;
    wav_props.s_R = wav_props.K_R - N + 2*max(Nprime,N0R); dual_wav_props.s_R = wav_props.s_R;
    wav_props.t_L = wav_props.s_L + max(wav_props.R-1,-dual_wav_props.L); 
    wav_props.t_R = wav_props.s_R + max(wav_props.R-1,-dual_wav_props.L);
    dual_wav_props.t_L = dual_wav_props.s_L + max(dual_wav_props.R-1,-wav_props.L);
    dual_wav_props.t_R = dual_wav_props.s_R + max(dual_wav_props.R-1,-wav_props.L);
    
    % Common DWT/dual DWT exceptions
    if wav_props.s_L + wav_props.s_R > dimphi1
        throw(MException('WL:bd_toomanylevels', 'Not enough room for all modified boundary functions at lowest resolution'));
    end
    
    % DWT exceptions
    if wav_props.t_L + Nprime > dimphi1
        throw(MException('WL:bd_toomanylevels', 'Expressions for left boundary functions need right boundary functions'));
    end
    if wav_props.t_R + Nprime > dimphi1
        throw(MException('WL:bd_toomanylevels', 'Expressions for right boundary functions need left boundary functions'));
    end
    
    % Dual DWT exceptions
    if dual_wav_props.t_L + Nprime > dimphi1
        throw(MException('WL:bd_toomanylevels', 'Expressions for dual left boundary functions need right boundary functions')); % TODO: Check these only in case of the dual transform
    end
    if dual_wav_props.t_R + Nprime > dimphi1
        throw(MException('WL:bd_toomanylevels', 'Expressions for dual right boundary functions need left boundary functions'));
    end
    
    % Mirror right hand side variables
    wav_props.A_R_pre_inv = fliplr(flipud(wav_props.A_R_pre_inv));
    dual_wav_props.A_R_pre_inv = fliplr(flipud(dual_wav_props.A_R_pre_inv));
    WR = fliplr(flipud(WR)); WRtilde = flipud(fliplr(WRtilde));
    if wav_props.L+wav_props.R==1
        for k = (size(WR,2)-(wav_props.K_R-N)):(-2):1
            WR(:, [k-1 k]) = WR(:, [k k-1]);
        end
        for k = (size(WRtilde,2)-(wav_props.K_R-N)):(-2):1
            WRtilde(:, [k-1 k]) = WRtilde(:, [k k-1]);
        end
    end
end

function [filtermatr,d]=compute_filter_matrix(wp, N, Nprime, length_signal)
    % in phibvals (left or right) a column will contain the values of a boundary function evaluated at [K-N,K+R-2]
    phibvals_left  = Gsegment(wp.phivals, (wp.L+1):(wp.R-1), (-N+1):(wp.R-1),           (2 - wp.K_L - wp.R):0, 0) * wp.CL; % TODO: Does not handle N~=Ntilde. 
    phibvals_right = Gsegment(wp.phivals, (wp.L+1):(wp.R-1), (-N+1+wp.L+wp.R):(wp.R-1), (2 - wp.K_R - wp.R):0, 0) * wp.CR;
    phibvals_right = phibvals_right( end:(-1):1, end:(-1):1); % Now contains values evaluated at M - [K-N+L+R,K+R-2].
    
    d = max( [-wp.L,N,wp.R]);
    modval = 2*d+1;
    filtermatr = zeros( modval, length_signal);
    
    % Set function values into filtermatr
    filtermatr( 1:(N+wp.R-1), 1:Nprime) = phibvals_left;
    for k=(Nprime+1):(length_signal- Nprime)
        filtermatr( mod(k - 1 + (( wp.L + 1):(wp.R - 1)), modval) + 1, k) = wp.phivals;
    end
    rowinds = mod( length_signal + ((wp.L - N + 2):0) - 1, modval) + 1;
    filtermatr( rowinds, (end - Nprime + 1):end ) = phibvals_right;
    
    filtermatr=luimpl_banded( filtermatr, d);
end

% The following functions assume that the bandwidth is d, and that the (i,j) entry is stored in A_{mod(i,2*d+1),j}
function A=luimpl_banded(A, d)
    n = size(A,2);
    for k=1:(n-1)
        currentrow = mod(k - 1, 2*d + 1) + 1;
        rows = mod( ((k+1):(min(k + d, n))) - 1, 2*d + 1) + 1;
        cols = (k+1):(min(k + d, n));
        A( rows, k)       = A(rows,k) /A(currentrow, k);
        A( rows, cols) = A(rows, cols) - A(rows, k)*A(currentrow, cols);
    end
end

% Following function also in bw_compute_left
function val=Gsegment(g0,supp,rowrange,colrange,symbolic)
    if symbolic 
        val = sym(zeros(length(rowrange),length(colrange)));
    else
        val = zeros(length(rowrange),length(colrange));
    end
    k=1;
    for col_ind = colrange
        actualinds =  supp + col_ind;
        [intersec,i1,i2] = intersect(rowrange,actualinds);
        val(i1,k) = g0( actualinds(i2) - actualinds(1) + 1 );
        k = k+1;
    end
end