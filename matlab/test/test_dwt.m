newmax = 0;
newmax=test_bd('spline4.4', 2, 255, 0, 4, 'filter', 'filter', newmax);
newmax=test_bd('spline2.6',    2, 255, 0, 2, 'filter', 'bd_pre', newmax);
newmax=test_bd('spline2.4',    2, 255, 0, 2, 'filter', 'bd_pre', newmax);
newmax=test_bd('spline4.4',    2, 255, 0, 4, 'filter', 'bd_pre', newmax);

test_dwt_different_sizes('spline4.4')

for k=0:1
    impl_strategy = 'filter';
    if k==0
        impl_strategy = 'lifting';
    end
    for dimx=127:129
        m = log2(64/(2*2));
        newmax=test_bd('spline53', m, dimx, 0, 2, impl_strategy, 'bd_pre', newmax);
        newmax=test_bd('pwl2',     m, dimx, 0, 2, impl_strategy, 'bd_pre', newmax);
    
        m = log2(64/(2*4));
        newmax=test_bd('cdf97',    m, dimx, 0, 4, impl_strategy, 'bd_pre', newmax);
    end

    for dimx=126:128
        m = log2(64/(2*4));
        newmax=test_bd('db4',      m, dimx, 1, 4, impl_strategy, 'bd_pre', newmax);
    end
end

compare_filter_lifting('db4', 'per')
compare_filter_lifting('spline53', 'symm')
compare_filter_lifting('cdf97', 'symm')

test_orthogonality('db2')
test_orthogonality('db4')
test_dwt_different_sizes('cdf97')
test_dwt_different_sizes('spline53')
test_dwt_different_sizes('pwl0')
test_dwt_different_sizes('pwl2')
test_dwt_different_sizes('haar')
test_dwt_different_sizes('spline4.4')
test_simple_dwt2()
test_haar()
newmax=test_bd_db_van(4, 'lifting', 'bd_pre', newmax);

newmax=test_bd('spline4.4', 2, 65, 0, 4, 'filter', 'bd_pre', newmax);
test_exceptions()

disp(sprintf('Highest deviation %0.5e', newmax))

function test_exceptions()
    disp('Checking that a lifting-based implementation for a wavelet that is not least dissimilar returns the expected error')
    x = rand(64,1);
    try
        x=wl_dwt_impl(x, 'spline4.4',  'bd_mode', 'bd', 'prefilter_mode', 'bd_pre', 'impl_strategy', 'lifting');
    catch exc
        assert( strcmpi(exc.identifier, 'WL:lifting_impossible') );
    end
    
    disp('Checking that the any mode defaults to filter-based implementation for the same wavelet')
    x=wl_dwt_impl(x, 'spline4.4',  'bd_mode', 'bd', 'prefilter_mode', 'bd_pre', 'impl_strategy', 'any');
end

function compare_filter_lifting(wave_name, bd_mode)
    disp(sprintf('Comparing lifting-based vs. filter-based: %s', wave_name))
    m = 4;
    x = rand(64,1);
    
    x1=wl_dwt_impl(x, wave_name, 'm', m, 'bd_mode', bd_mode, 'impl_strategy', 'filter');
    x2=wl_dwt_impl(x, wave_name, 'm', m, 'bd_mode', bd_mode, 'impl_strategy', 'lifting');
    diff = max(abs(x1-x2));
    assert(diff ~= 0 && diff < 1E-11)
    
    x1=wl_idwt_impl(x1, wave_name, 'm', m, 'bd_mode', bd_mode, 'impl_strategy', 'filter');
    x2=wl_idwt_impl(x2, wave_name, 'm', m, 'bd_mode', bd_mode, 'impl_strategy', 'lifting');
    diff = max(abs(x1-x2));
    assert(diff ~= 0 && diff < 1E-11)
    
    diff = max(abs(x-x1));
    assert(diff ~= 0 && diff < 1E-11)
    
    diff = max(abs(x-x2));
    assert(diff ~= 0 && diff < 1E-11)
end

function test_haar()
    disp('Testing Haar vs. filter-based Haar')
    m = 4;
    
    sd = rand(64,1);
    sd1 = wl_dwt_impl(sd, 'haar', 'm', m);
    
    [f, prefilter, offset_L, offset_R]=wl_find_kernel_from_filters([1 1]/sqrt(2), 1, [1 1]/sqrt(2), 1, 64, 1, 'm', m, 'impl_strategy', 'lifting');
    sd2=wl_dwt1_impl_internal(sd, f, prefilter, [offset_L offset_R], 'm', m, 'bd_mode', 'per');
    
    diff = max(abs(sd2-sd1));
    assert(diff ~= 0 && diff < 1E-13)
end

function test_dwt_different_sizes(wave_name)
    disp('Testing the DWT on different input sizes')
    m = 4;

    disp(sprintf('Testing 2D with one channel: %s', wave_name))
    img = rand(64);
    img2 = img;
    img2 = wl_dwt_impl(img2, wave_name, 'm', m, 'dims', 2);
    img2 = wl_idwt_impl(img2, wave_name, 'm', m, 'dims', 2);
    diff = max(max(abs(img2-img)));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 2D with three channels: %s', wave_name))
    img = rand(64, 64, 3);
    img2 = img;
    img2 = wl_dwt_impl(img2, wave_name, 'm', m);
    img2 = wl_idwt_impl(img2, wave_name, 'm', m);
    diff = max(max(max(abs(img2-img))));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 1D with one channel: %s', wave_name))
    sd = rand(64,1);
    sd2 = sd;
    sd2 = wl_dwt_impl(sd2, wave_name, 'm', m);
    sd2 = wl_idwt_impl(sd2, wave_name, 'm', m);
    diff = max(abs(sd2-sd));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 1D with two channels: %s', wave_name))
    sd = rand(64,2);
    sd2 = sd;
    sd2 = wl_dwt_impl(sd2, wave_name, 'm', m);
    sd2 = wl_idwt_impl(sd2, wave_name, 'm', m);
    diff = max(max(abs(sd2-sd)));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 3D with one channel: %s', wave_name))
    sd = rand(64,64,64);
    sd2 = sd;
    sd2 = wl_dwt_impl(sd2, wave_name, 'm', m, 'dims', 3);
    sd2 = wl_idwt_impl(sd2, wave_name, 'm', m, 'dims', 3);
    diff = max(max(max(abs(sd2-sd))));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 3D with three channels: %s', wave_name))
    sd = rand(64,64,64,3);
    sd2 = sd;
    sd2 = wl_dwt_impl(sd2, wave_name, 'm', m);
    sd2 = wl_idwt_impl(sd2, wave_name, 'm', m);
    diff = max(max(max(max(abs(sd2-sd)))));
    assert(diff ~= 0 && diff < 1E-13)
end

function test_orthogonality(wave_name)
    disp('Testing orthonormal wavelets:')
    x0 = rand(32,1);
    
    disp(sprintf('Testing that the IDWT inverts the DWT: %s', wave_name))
    x = x0;
    x = wl_dwt_impl(x, wave_name, 'm', 2);
    x = wl_idwt_impl(x, wave_name, 'm', 2);
    diff = max(abs(x-x0));
    assert(diff ~= 0 && diff < 1E-13)

    disp(sprintf('See that the wavelet transform equals the dual wavelet transform: %s', wave_name))
    x = x0;
    x = wl_dwt_impl(x, wave_name, 'm', 2, 'dual', 1);
    x0 = wl_dwt_impl(x0, wave_name, 'm', 2, 'dual', 0);
    diff = max(abs(x-x0));
    assert(diff ~= 0 && diff < 1E-13)

    disp(sprintf('Apply the transpose, to see that the transpose equals the inverse: %s', wave_name))
    x = x0;
    x = wl_dwt_impl(x, wave_name, 'm', 2, 'transpose', 1);
    x = wl_dwt_impl(x, wave_name, 'm', 2, 'transpose', 0);
    diff = max(abs(x-x0));
    assert(diff ~= 0 && diff < 1E-13)
end

function test_simple_dwt2()
    disp('Testing simple DWT2')
    img = rand(32, 32, 3);
    img2 = img;
% Begin simple_dwt2
    f = @(x, bd_mode) wl_dwt_impl(x, 'cdf97', 'm', 4, 'bd_mode', bd_mode, 'dims', 1);
    img = tensor2_impl(img, f, f, 'symm');
% End simple_dwt2
% Begin simple_idwt2
    invf = @(x, bd_mode) wl_idwt_impl(x, 'cdf97', 'm', 4, 'bd_mode', bd_mode, 'dims', 1);
    img = tensor2_impl(img, invf, invf, 'symm');
% End simple_idwt2
    diff = max(max(max(abs(img2-img))));
    assert(diff ~= 0 && diff < 1E-13)
end

function newmax=test_bd(wave_name, m, dimx, L_p_R, N, impl_strategy, prefilter_mode, newmax)
    disp(sprintf('Testing bd %s with dimx=%i. %s-based. prefilter mode %s', wave_name, dimx, impl_strategy, prefilter_mode));
    
    res2 = dimx + L_p_R - 1;
    res3 = mod(res2, 2^m);
    if res3 == 0 
        toadd = 0;
    else
        toadd = 2^m - res3;
    end
    dimphi0 =  2^(-m)*dimx + (1-2^(-m))*(1 - toadd - L_p_R);
    
    for s=0:(N-1)
        res = ((1:dimx).^s)';
        x=wl_dwt_impl(res, wave_name,  'm', m, 'bd_mode', 'bd', 'prefilter_mode', prefilter_mode, 'impl_strategy', impl_strategy);
        maxval = max(abs(x((dimphi0+1):end))); newmax=max([newmax,maxval]);
        disp(sprintf('highest detail: %0.5e', maxval));
        x=wl_idwt_impl(x,  wave_name, 'm', m, 'bd_mode', 'bd', 'prefilter_mode', prefilter_mode, 'impl_strategy', impl_strategy);
        diff = max(abs(res-x)); newmax=max([newmax,diff]);
        disp(sprintf('highest diff: %0.5e', diff));
    end
end


function newmax=test_bd_db_van(N, impl_strategy, prefilter_mode, newmax)
    m=2;
    for k=2:N
        newmax=test_bd(sprintf('db%i',k), m, 64, 1, k, impl_strategy, prefilter_mode, newmax);
    end
end

% x = double(imread('images/lena.png'));
% x=dwt_impl(x, 'cdf97', 'm', 1, 'bd_mode', 'bd', 'prefilter_mode', 'bd_pre'); 
% x=idwt_impl(x, 'cdf97', 'm', 1, 'bd_mode', 'bd', 'prefilter_mode', 'bd_pre');
% imshow(uint8(x))