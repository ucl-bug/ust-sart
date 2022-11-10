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
- `temperature`: temperature of the water during data acquisition (degrees C).

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

## Upsampling

For fast computation, it is common to start the image reconstruction on a coarse grid, and to then upsample the computational grid as the iterations progress. In `ust-sart`, the array `ups` is defined which contains the upsample factors for each iteration, relative to the starting grid size. For example:
```
ups        = [1, 1, 2, 2, 4, 4]; % upsampling factors for each iteration
Nit        = length(ups);  % number of iterations
dx0        = 8e-3;         % step size for iteration 1 [m]
``` 
In this case, six iterations would be performed, starting with a step size of 8 mm, progressing to a step size of 4 mm, and finishing with a step size of 2 mm. For every iteration, the physical length of the computational grid remains the same, and the sound speed map is interpolated to match the required grid step size.

## Getting Started

Clone this repository in a terminal: `git clone https://github.com/ucl-bug/ust-sart.git`
Open Matlab2022a or more recent, and add the folder to the path:
```
addpath(genpath('<pathname>\ust-sart'));
```
This code uses the functions `waterSoundSpeed` and `getWin` from the [k-Wave toolbox](https://github.com/ucl-bug/k-wave).
Change directories to the example folder, and run the example script `sart_example.m`.

```
close all
clearvars

% Make sure ust-sart is on the path
addpath(genpath('..\')); 
load('example_element_positions');
load('example_delta_tof');

% remove lower triangle TOF data (it's reciprocal) to save computation time
delta_tof(logical(tril(ones(size(delta_tof)), -1))) = NaN;

% Initialise sart object
temperature = 20;   % water temperature [degC]
sart        = SartExperiment(element_positions, delta_tof, temperature);

% Plot setup and data
sart.plotSetup;

% Perform reconstruction
ups        = [1, 1, 2, 4]; % upsampling factors for each iteration
Nit        = length(ups);  % number of iterations
dx0        = 4e-3;         % step size for iteration 1 [m]
init_c_val = sart.c_water; % sound speed value for homogeneous initial estimate [m/s]
hamming    = 0;            % boolean controlling whether hamming window is used
recon_d    = 0.135;        % diameter of reconstruction circle [m]
sart.reconstructSart(init_c_val, recon_d, Nit, dx0, ups, hamming=hamming);

% plot final estimate
figure; imagesc(squeeze(sart.estimates(:,:,end))); axis image;
```
![setup_example](https://github.com/ucl-bug/ust-sart/blob/main/setup_example.png)
![recon_example](https://github.com/ucl-bug/ust-sart/blob/main/recon_example.png)