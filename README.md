# ust-sart
Code for reconstructing 2D sound speed maps from time-of-flight data using straight-ray-SART.

This code is intended to be used to reconstruct ultrasound tomography (UST) data acquired from 2D ring array transducer hardware, such as the [open-UST](https://github.com/morganjroberts/open-UST) system.

**Note:** This code does not handle absorption or density, and only inverts for sound speed.
**Note:** This code does not provide functions for time-of-flight picking.

## Background

- This inverse problem cannot be solved by direct matrix inversion because the coefficient matrix is too large to be stored in memory.
- Alternatively, iterative methods should be used, such as Kaczmarz's method of projections.
- This code is an implementation of the straight ray tracing technique described by Kak and Slaney in Chapter 7 of `A. C. Kak and Malcolm Slaney, Principles of Computerized Tomographic Imaging, IEEE Press, 1988` (section 7.4, or p.11 in the pdf), available [here](https://www.slaney.org/pct/pct-toc.html), or in this repository.

## Required input data
For a transducer ring array with N elements:

- `delta_tof`: this is a N x N matrix indexed as `(tdx, rdx)`. Each element contains the time-of-flight difference in seconds between the phantom UST data and the watershot UST data: `delta_tof(tdx, rdx) = tof_phantom(tdx, rdx) - tof_water(tdx, rdx)`.
- `element_positions`: this is a N x 2 matrix containing the cartesian coordinates of each transducer element. The ring array should be centred on `(0, 0)`.

## Example
