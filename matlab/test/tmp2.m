clear all;

fname = 'lena.png';

im = double(rgb2gray(imread(fname)));
size(im)

m = 2;
%prefilter_mode = 'bd_pre';
wname = 'cdf97';
dims = 2;

bd_mode = 'symm';
prefilter_mode = 'none';
dwt_sym = wl_dwt_impl(im, wname, 'm', m, 'bd_mode', bd_mode, 'prefilter_mode', prefilter_mode, 'dims', dims);

bd_mode = 'bd';
prefilter_mode = 'bd_pre';
dwt_bd =  wl_dwt_impl(im, wname, 'm', m, 'bd_mode', bd_mode, 'prefilter_mode', prefilter_mode, 'dims', dims);
inv    = wl_idwt_impl(dwt_bd, wname, 'm', m, 'bd_mode', bd_mode, 'prefilter_mode', prefilter_mode, 'dims', dims);
max(max(abs(inv-im)))

N = size(im,1);
im_sym = dwt_sym(N/2+1:end, N/2+1:end);
im_bd  = dwt_bd(N/2+1:end, N/2+1:end);

mi = min(im_sym(:))
ma = max(im_sym(:))

mibd = min(im_bd(:))
mabd = max(im_bd(:))

% To get same colorbar on both images
im_bd(im_bd<mi) = mi;
im_bd(im_bd>ma) = ma;

figure();
subplot(121); imagesc(abs(im_sym)); colormap('gray'); colorbar(); title('symmetric');
subplot(122); imagesc(abs(im_bd)); colormap('gray'); colorbar(); title('VMP')


