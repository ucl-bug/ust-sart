# ust-sart
Code for reconstructing 2D sound speed maps from time-of-flight data using straight-ray-SART.

This code is intended to be used to reconstruct ultrasound tomography (UST) data acquired from 2D ring array transducer hardware, such as the [open-UST](https://github.com/morganjroberts/open-UST) system.

**Note:** This code does not handle absorption or density, and only inverts for sound speed.

**Note:** This code does not provide functions for time-of-flight picking.

## Background

- This inverse problem cannot be solved by direct matrix inversion because the coefficient matrix is too large to be stored in memory.
- Alternatively, iterative methods should be used, such as Kaczmarz's method of projections.
- This code is an implementation of the straight ray tracing technique described by Kak and Slaney in Chapter 7 of `A. C. Kak and Malcolm Slaney, Principles of Computerized Tomographic Imaging, IEEE Press, 1988` (section 7.4, or p.11 in the pdf), available [here](https://www.slaney.org/pct/pct-toc.html), or in this repository.

## Required Input Data
For a transducer ring array with N elements:

- `delta_tof`: this is a N x N matrix. Each `(tdx, rdx)` value contains the time-of-flight difference in seconds between the phantom UST data and the watershot UST data, for the ray joining the `tdx` and `rdx` transducer elements: `delta_tof(tdx, rdx) = tof_phantom(tdx, rdx) - tof_water(tdx, rdx)`. Values of `NaN` should be used for missing data, not `0`.
- `element_positions`: this is a N x 2 matrix containing the cartesian coordinates of each transducer element. The ring array should be centred on `(0, 0)`.

**Note:** Reciprocity means that a ray joining `tx = 1` with `rx = 139` should have the same time-of-flight data as a ray joining `tx = 139` with `rx = 1`. Therefore, for faster computation the lower triangular part of `delta_tof` should be set to `NaN` so that the code ignores these rays.

## Slowness and Sound Speed
For every iteration, the sound speed estimate `c_est` is first converted to relative slowness:
```
s_est = ((1 ./ c_est) - (1 ./ c_water));
```
The slowness difference relative to water `s_est` is then used to solve the effective linear alebgra problem 
```
delta_tof = d * s_est;
```
using Kaczmarz's method of projections, since the distance matrix `d` is too large to store in memory.
Finally, the new slowness estimate `new_s_est` is converted back to sound speed:
```
new_c_est = c_water ./ ((new_s_est .* c_water) + 1);
```
This process repeats for each iteration.

## Getting Started

Clone this repository in a terminal: `git clone https://github.com/ucl-bug/ust-sart.git`
Open Matlab2022a or more recent, and add the folder to tha path:
```
addpath('<pathname>/ust-sart');
```

Then, run the example script `sart_example.m`.

```
close all
clearvars

load('256_element_UST_array_predicted_positions.mat', 'element_positions');
load('delta_tof_phantom.mat', 'delta_tof');

% remove lower triangle data (it's reciprocal) to save computation time
delta_tof(logical(tril(ones(size(delta_tof)), -1))) = NaN;

% step size for iteration 1
dx0 = 4e-3; 
% upsampling factor for each iteration (relative to original step size)
upsample_factors = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4, 4, 4]; 
Nit = length(upsample_factors);

sart = SartExperiment;
sart.deriveParams(dx0);
sart.loadDetectorPositions(element_positions)
sart.loadTOFdata(delta_tof, 20)
sart.plotSetup;

init_est = 1480 * ones(sart.Nx);
hamming = 0; % boolean (hamming window used to backpropagate errors for each ray?)
sart.reconstructSart(init_est, Nit, upsample_factors, hamming);
```

![recon_example](https://github.com/ucl-bug/ust-sart/blob/main/recon_example.png)