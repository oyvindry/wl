# wl - General purpose wavelet library

This repository contains the source code for the wl wavelet library. 
It contains open source implementations of the discrete wavelet transform (DWT).
In particular we highlight the advantages of the implementation.

* Filter-based and lifting-based implementations (which one to use can be specified by the caller).
* Support for wavelets preserving vanishing moments on the interval (not only for Spline wavelets and for orthonormal wavelets). 
  To the best of our knowledge this is the only openly available software implementation which computes boundary wavelets on the fly, and which supports general input sizes.

## Installation

### MATLAB
Add the folder `MATLAB` to your MATLAB path. To do this automatically on startup, add the folder to your 
[startup.m file](https://ch.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html). 

## User guide to the discrete wavelet transform implementation

The software has two main functions `wl_dwt_impl` and `wl_idwt_impl` which implement
a DWT and a IDWT, respectively. These functions can be called with a variety of
arguments, specifying the type of wavelet transform and its properties. To
compute an `m`-level DWT on a vector `x`, using a wavelet with name `wname`,
boundary handling mode `bd_mode`, and prefiltering mode `prefilter_mode`,
simply write
```
y = wl_dwt_impl(x, wname, 'm', m, 'bd_mode', bd_mode, 'prefilter_mode', prefilter_mode);
```
Here all but the two first arguments are optional. If they are not provided the implementaion will set maningful default values for them, as described in the documentation for the function `wl_setopts`. 
Here `wname` will be one of the following supported wavelets:
*  `'cdf97'` : CDF 9/7 wavelet with $N = \tilde{N} = 4$ vanishing moments.
*  `'spline53'` : Spline 5/3 wavelet with $N = \tilde{N} = 2$ vanishing moments.
*  `'splineX.X'` : Spline wavelet with `X` number of vanishing moments for the wavelet and dual wavelet,
*  `'pwl0'` : Piecewise linear wavelets with 0 vanishing moments,
*  `'pwl2'` : Piecewise linear wavelets with 2 vanishing moments,
*  `'Haar'` : The Haar wavelet,
*  `'dbX'` : Daubechies orthonormal wavelet with `X` vanishing moments,
*  `'symX'` : Symlets. A close to symmetric, orthonormal (Daubechies) wavelet with `X` vanishing moments.

Likewise `bd_mode` can take the values
* `'per'` : Periodic extension,
* `'symm'` : Symmetric extension (not for orthonormal wavelets),
* `'none'` : No boundary handling, in the sense that the input is zero-padded. 
* `'bd'` : Boundary handling preserving vanishing moments on the interval.

Note that, when the `'bd'` mode is invoked, computations in the `'none'` mode are also performed by the software, in order to compute the tail handling components. 
`prefilter_mode` can take the values
* `'none'` : No prefiltering (default),
* `'bd_pre'` : Boundary wavelets with preconditioning as described in XXXX (Only available in MATLAB version) 

Not all combinations of these arguments make sense. For instance it is not possible to apply a symmetric boundary extension to an orthonormal
wavelet. In such cases the functions halt, issuing an error. 

`wl_dwt_impl` also accepts the following arguments.
* `dims` : The number of dimensions to apply the DWT to. If the input is two-dimensional, this enables the caller to specify whether a two dimensional DWT should be applied, or a one dimenionsal DWT vectorized on the second axis. 
* `dual` :  Whether the dual transform should be applied or not.
* `transpose` : Whether the transform shoudl be transposed.
* `data_layout` : How data should be assembled. Possible values are `resolution` (lowest resolution first, default), and `time` (sort according to time). 
    
    
    
## Internal functions and efficient computations 
The `wl_dwt_impl` and `wl_idwt_impl` functions compute the
filter coefficients and tail handling components on the fly each time the
functions are invoked. This allows the software to compute different boundary
functions for different input sizes and makes the functions user friendly. At
the same time, this increases the computational time as each call recomputes coefficients. To allow for using precomputed boundary functions, the software has the functions `wl_dwt1_impl_internal` and
`wl_idwt1_impl_internal` which will (given the right input) compute a one dimensional DWT and IDWT, respectively, using precomputed boundary functions.  Similar functions exist for two and three
dimensions. 

To use these internal functions we need to first compute the filter coefficients 
manually. The complete code can be as follows.
```
[f, prefilter, offset_L, offset_R]=wl_find_kernel(wave_name, length_signal, 1)
save('wl_kernel.mat', 'f', 'prefilter', 'offset_L', 'offset_R');
...
load('wl_kernel.mat');
x = wl_dwt1_impl_internal(x, f, prefilter, [offset_L offset_R]);
```
Let us go through the different pieces of this code. 
* The `wl_find_kernel` functions return two function handles:
  1. `f(x, bd_mode)` performs a one level in-place DWT on `x` with boundary mode `bd_mode`, 
  2. `prefilter(x, forward)` filters `x`, where `forward` is either 0 (postfiltering) or 1 (prefiltering). 
* The `offset` parameter is only necessary when using boundary wavelets. It tells `wl_dwt1_impl_internal` that there are $N-K_L$ and $N-K_R$ extra wavelets at each boundary.

Since the internal functions avoid precomputation, their execution times 
should be comparable with those of convolution, since the DWT/IDWT can be 
expressed in terms of this (and possibly lifting).

Following the code above, it is possible to experiment with custom made kernels - simply implement your own kernel function (taking the same arguments as the ones defined in the software), and use it as input to `wl_dwt1_impl_internal`.

## Citing
When referencing the wl library, please consider to cite the following two references

```
@article{antun2021unification,
  title={On the Unification of Schemes and Software for Wavelets on the Interval},
  author={Antun, Vegard and Ryan, {\O}yvind},
  journal={Acta Applicandae Mathematicae},
  volume={173},
  number={1},
  pages={1--25},
  year={2021},
  publisher={Springer}
}

@book{ryan2019linear,
  title={Linear Algebra, Signal Processing, and Wavelets--a Unified Approach},
  author={Ryan, {\O}yvind and Ryan and Peters},
  year={2019},
  publisher={Springer}
}
```

