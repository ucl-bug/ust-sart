function delta_d_ijm = interpolateRaySection(obj, m_xvec, m_yvec, mdx)
%INTERPOLATERAYSECTION Calculate the contribution of neighbouring pixels to a ray section
%
% DESCRIPTION:
%     interpolateRaySection locates a small section of a ray joining a
%     tdx/rdx. It locates the neighbouring pixels in the sound speed map
%     and calculates how much each pixel contributes to the sound
%     propagation through this small ray section, using bilinear
%     interpolation. The basis functions obtained from bilinear elements
%     are pyramid shaped, each with a support extending over a square
%     region the size of four pixels. The bilinear interpolation allows a
%     continuous form of the underlying sound speed map to be regenerated
%     for comnputation. This function is designed to be called iteratively
%     across all ray sub-sections for mdx = 1, 2, ... M.
%
% USAGE:
%     delta_d_ijm = interpolateRaySection(obj, m_xvec, m_yvec, mdx)
%
% INPUTS: 
%     obj       - object, instance of the SartExperiment class
%     m_xvec    - [numeric] M x 1 array, ray sections x-coordinates [m]
%     m_yvec    - [numeric] M x 1 array, ray sections y-coordinates [m]
%     mdx       - [numeric] integer index selecting the current ray section
%
% OUTPUTS:
%     delta_ijm - [numeric] square matrix the same size as the current
%                 sound speed map, containing the pixel-wise contribution
%                 to the current ray section. Has values between 0-1.
%
% REFERENCES:
%     A. C. Kak and Malcolm Slaney, Principles of Computerized Tomographic
%     Imaging, IEEE Press, 1988 (section 7.4), available at
%     https://www.slaney.org/pct/pct-toc.html
%
% ABOUT:
%     author    - Morgan Roberts
%     date      - 9th November 2022


% find coordinates of neighbouring grid points (checking for points within a spatial step size)
mx = m_xvec(mdx);
my = m_yvec(mdx);
Ix = find(abs(obj.grid_x - mx) < obj.dx);
Iy = find(abs(obj.grid_x - my) < obj.dx);

% catch for case where ray section located exactly on a grid
% point
if length(Ix) == 1
    Ix = [Ix, Ix + 1];
end
if length(Iy) == 1
    Iy = [Iy, Iy + 1];
end

% error checking
if isempty(Ix)
    error('Ix is empty, could not find nearest grid point')
end

% extract coordinates of the neighbouring grid points
x1 = obj.grid_x(Ix(1));
x2 = obj.grid_x(Ix(2));
y1 = obj.grid_x(Iy(1));
y2 = obj.grid_x(Iy(2));

% find coefficients of bilinear interpolation b
% https://en.wikipedia.org/wiki/Bilinear_interpolation
mat = [1, x1, y1, x1 * y1;
       1, x1, y2, x1 * y2;
       1, x2, y1, x2 * y1;
       1, x2, y2, x2 * y2];
vec = [1; mx; my; mx * my];
b   = inv(mat)' * vec;

% check for negative weights
if sum(b < -eps*1e3) > 0
    error('negative weights in d_vec');
end

% locate weights within domain matrix, vectorise and cumulative sum
% with d_ijm vector
temp_mat = zeros(obj.Nx);
temp_mat(Ix(1), Iy(1)) = b(1);
temp_mat(Ix(1), Iy(2)) = b(2);
temp_mat(Ix(2), Iy(1)) = b(3);
temp_mat(Ix(2), Iy(2)) = b(4); 

% adjust first/last section so total physical length of ray preserved
M = length(m_xvec);
if mdx == 1 || mdx == M
    temp_mat = temp_mat * 0.5;
end

delta_d_ijm = temp_mat;

end