function [f, prefilter, offset_L, offset_R]=wl_find_kernel(wave_name, length_signal, forward, varargin)
    % Computes the properties of a wavelet with the given name. What properties 
    % are computed depend on the bd_mode parameter, m, and length_signal.
    %
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
    % length_signal: Length of the input signal. Default: 0.
    % forward: Whether the forward transform should be used
    %
    % This function also accepts a number of named, optional parameters. These are parsed by the function wl_setopts(). 
    % The documentation of this function also contains the full documentation of these optional parameters.
    
    args = {'wave_name', wave_name, varargin{:}};
    opts = wl_setopts(args{:});
    
    if strcmpi(wave_name,'Haar')
        prefilter = @(x, forward) x;
        if opts.transpose
            forward = ~forward;
        end
        if forward
            f = @dwt_kernel_haar;
        else
            f = @idwt_kernel_haar;
        end
        offset_L = 0; offset_R = 0;
        return;
    end
    
    N = 0;
    Ntilde = 0;
    if strcmpi(wave_name(1:2), 'db')
        N = str2double(wave_name(3:end)); Ntilde = N;
        [g0,h0]=compute_ortho(N);
    elseif strcmpi(wave_name(1:3), 'sym')
        N = str2double(wave_name(4:end)); Ntilde = N;
        [g0,h0]=compute_sym(N);
    elseif strcmpi(wave_name, 'pwl0')
        g0 = [1/2 1 1/2]/sqrt(2);
        h0 = sqrt(2);
    elseif strcmpi(wave_name, 'pwl2')
        N = 2; Ntilde = 2;
        g0 = [1/2 1 1/2]/sqrt(2);
        h0 = [-1/8 1/4 3/4 1/4 -1/8]*sqrt(2);
    elseif strcmpi(wave_name, 'spline53')
        N = 2; Ntilde = 2;
        g0 = [1/4 1/2 1/4];
        h0 = [-1/4 1/2 3/2 1/2 -1/4];
    elseif strcmpi(wave_name, 'cdf97')
        N = 4; Ntilde = 4;
        [g0,h0]=compute_97();
    elseif strcmpi(wave_name(1:6), 'spline')
        N = str2double(wave_name(7));
        Ntilde = str2double(wave_name(9));
        [g0,h0]=compute_spline_filters(N, Ntilde);
    else
        throw(MException('WL:wavelet_not_supported', 'Wavelet not supported'));
    end
    
    [f, prefilter, offset_L, offset_R]=wl_find_kernel_from_filters(g0, N, h0, Ntilde, length_signal, forward, args{:});
end


function x=dwt_kernel_haar(x, bd_mode)
    x = x/sqrt(2);
    N = size(x, 1);
    if mod(N,2) == 1
        x(1:2, :) = [x(1, :) + x(2, :) - x(N, :);x(1, :) - x(2, :) - x(N, :)];
        x(N, :) = 2*x(N, :);
    else
        x(1:2, :) = [x(1, :) + x(2, :); x(1, :) - x(2, :)];
    end
    for k = 3:2:(N-1)
        x(k:(k+1), :) = [x(k, :) + x(k+1, :); x(k, :) - x(k+1, :)];
    end
end
% End dwt_kernel_haar

function x=idwt_kernel_haar(x, bd_mode)
    x = x/sqrt(2);
    N = size(x, 1);
    if mod(N,2) == 1
        x(1:2, :) = [x(1, :) + x(2, :) + x(N, :); x(1, :) - x(2, :)];
    else
        x(1:2, :) = [x(1, :) + x(2, :); x(1, :) - x(2, :)];
    end
    for k = 3:2:(N-1)
        x(k:(k+1), :) = [x(k, :) + x(k+1, :); x(k, :) - x(k+1, :)];
    end  
end
% End idwt_kernel_haar

function [g0,h0]=compute_spline_filters(N, Ntilde)
  Navg=(N+Ntilde)/2;
  vals=compute_QN(Navg);
  
  h0 = 1;
  for k=1:(Ntilde/2)
    h0=conv(h0,[1/4 1/2 1/4]);
  end
  h0 = conv(h0, vals);
  
  g0=1;
  for k=1:(N/2)
    g0=conv(g0,[1/4 1/2 1/4]);
  end
  
  %h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2));
  %g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
end

function [g0,h0]=compute_97()
  N=4;
  vals=compute_QN(N);
   
  rts=roots(vals)';
  rts1=rts(find(abs(imag(rts))>0.001)); % imaginary roots
  rts2=rts(find(abs(imag(rts))<0.001)); % real roots
  
  h0=1;
  for rt=rts1
    h0=conv(h0,[-rt 1]);
  end
  for k=1:(N/2)
    h0=conv(h0,[1/4 1/2 1/4]);
  end
  h0=h0*vals(1);
  
  g0=1;
  for rt=rts2
    g0=conv(g0,[-rt 1]);
  end
  for k=1:(N/2)
    g0=conv(g0,[1/4 1/2 1/4]);
  end
  
  
  h0=real(h0);
  g0=real(g0);  
  x = sqrt(2)/abs(sum(h0));
  g0=g0/x;
  h0=h0*x;
  g0 = (g0+flip(g0))/2;
  h0 = (h0+flip(h0))/2;
  
  %h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2));
  %g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
end

function [g0,h0]=compute_ortho(N)
    % Computes the wavelet coefficients of the orthonormal Daubechies wavelet
    % N vanishing moments and with minimum phase   
    vals=compute_QN(N);
    rts=roots(vals)';
    rts1=rts(find(abs(rts)>1));

    g0=1;
    for rt=rts1
        g0=conv(g0,[-rt 1]);
    end
    g0 = real(g0);
    K=sqrt(vals(1)*(-1)^(length(rts1))/abs(prod(rts1)));
    g0=K*g0;
    for k=1:N
        g0=conv(g0,[1/2 1/2]);
    end
    
    % Ensuring integral is positive - This part of the code requiere some more
    % testing
    if (sum(g0) < 0)
        g0 = -g0;
    end
    h0=fliplr(g0);
    %g1=h0.*(-1).^(0:(length(g0)-1)); % Should be multiplied by -1? 
    %h1=g0.*(-1).^(1:(length(g0))); % Should be multiplied by -1? 
    %%h1=fliplr(g1);
end

function [g0,h0]=compute_sym(N)
    % Computes the wavelet coefficients of the orthonormal wavelet with N
    % vanishing moments and close to linear phase. This makes the wavelet
    % almost symmetric. These wavelets are called 'symlets'
    %
    % This function relies on matlabs wavelet coefficients. In the next version
    % this will be changed 
    
    currDWTmode = dwtmode('status', 'nodisp');
    dwtmode('per','nodisp');
    nu = 7;
    n = 2^nu;
    x = zeros([1,n]);
    x(ceil(N/2)) = 1;
    
    S = [2^(nu-1); 2^(nu-1); n]; % compute the S given by wavedec
    wave_name = sprintf('sym%d', N);
    
    y = waverec(x, S, wave_name);
    if (mod(N,2) == 1) % is odd
        g0 = y(1:2*N);
    else % is even 
        g0 = [y(end), y(1:2*N-1)];
    end
    
    h0=fliplr(g0);
    %g1=h0.*(-1).^(0:(length(g0)-1)); 
    %h1=fliplr(g1);
    
    dwtmode(currDWTmode, 'nodisp');
end

function vals=compute_QN(N)
    % Compute the coefficients in Q^(N)((1-cos(w))/2)
    k=0:(N-1);
    QN = 2*factorial(N+k-1)./(factorial(k).*factorial(N-1));
    vals=zeros(1,2*N-1);
    vals=QN(1);
    start=1;
    for k=2:N
        start=conv(start,[-1/4 1/2 -1/4]);
        vals=[0 vals 0]+QN(k)*start;
    end
end