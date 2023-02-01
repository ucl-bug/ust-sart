function [d_ijm, t_ij] = calculatePixelWeights(obj, m_xvec, m_yvec)
%CALCULATEPIXELWEIGHTS Calculate the contribution of total pixels to a ray
%
% DESCRIPTION:
%     calculatePixelWeights takes a set of M ray-section coordinates,
%     locates the neighbouring pixels in the sound speed map and calculates
%     how much each pixel contributes to the sound propagation through each
%     ray section, using bilinear interpolation. The basis functions
%     obtained from bilinear elements are pyramid shaped, each with a
%     support extending over a square region the size of four pixels. The
%     bilinear interpolation allows a continuous form of the underlying
%     sound speed map to be regenerated for comnputation. The result d_ijm
%     contains the cumulative pixel weightings for all the ray sections. A
%     hamming-weighted version of the result is also returned, which can be
%     used to favour sound speed corrections nearer to the centre of the
%     grid. This function should be called once for each ray.
%
% USAGE:
%     [d_ijm, t_ij] = calculatePixelWeights(obj, m_xvec, m_yvec)
%
% INPUTS: 
%     obj       - object, instance of the SartExperiment class
%     m_xvec    - [numeric] M x 1 array, ray sections x-coordinates [m]
%     m_yvec    - [numeric] M x 1 array, ray sections y-coordinates [m]
%
% OUTPUTS:
%     d_ijm     - [numeric] square matrix the same size as the current
%                 sound speed map, containing the pixel-wise contribution
%                 to the ray section.
%     t_ij      - [numeric] square matrix the same size as the current
%                 sound speed map, containing the hamming-weighted
%                 pixel-wise contribution to the ray section.
%
% REFERENCES:
%     A. C. Kak and Malcolm Slaney, Principles of Computerized Tomographic
%     Imaging, IEEE Press, 1988 (section 7.4), available at
%     https://www.slaney.org/pct/pct-toc.html
%
% ABOUT:
%     author    - Morgan Roberts
%     date      - 1st February 2023

M     = length(m_xvec);
d_ijm = zeros((obj.Nx .^2), 1, 'double'); 
t_ij  = 0;

% Locate the grid points neighbouring the ray section centroids
px         = abs(m_xvec - obj.grid_x') < obj.dx;
py         = abs(m_yvec - obj.grid_x') < obj.dx;
px         = diff(px,1,1);
py         = diff(py,1,1);
px         = px > 0;
py         = py > 0;
[i1,~]     = find(px);
[j1,~]     = find(py);
i1         = i1 + 1;
j1         = j1 + 1;
i2         = i1 + 1;
j2         = j1 + 1; % i2, j2 shouldn't exceed Nx since padding is used
if length(px) ~= length(py)
    warning('When locating neighbouring grid points, px and py different lengths');
end

% Coordinates of the neighbouring grid points, for each ray section
x1 = obj.grid_x(i1);
x2 = obj.grid_x(i2);
y1 = obj.grid_x(j1);
y2 = obj.grid_x(j2);

% calculate pixel weights with bilinear interpolation, for each ray section
denominator = (x2 - x1) .* (y2 - y1);
w11 = ((x2 - m_xvec) .* (y2 - m_yvec)) ./ denominator;
w12 = ((x2 - m_xvec) .* (m_yvec - y1)) ./ denominator;
w21 = ((m_xvec - x1) .* (y2 - m_yvec)) ./ denominator;
w22 = ((m_xvec - x1) .* (m_yvec - y1)) ./ denominator;
if sum( (w11 + w12 + w21 + w22)  - 1 > 1e-6 ) > 0
    error('Pixel weights do not sum to 1 for current ray');
end

% Adjust first/last section so total physical length of ray preserved
w11([1, end]) = w11([1, end]) * 0.5;
w12([1, end]) = w12([1, end]) * 0.5;
w21([1, end]) = w21([1, end]) * 0.5;
w22([1, end]) = w22([1, end]) * 0.5;

% Calculate the linear index of the pixel weights for each ray section
sz = [obj.Nx, obj.Nx];
k11 = sub2ind(sz, i1, j1);
k12 = sub2ind(sz, i1, j2);
k21 = sub2ind(sz, i2, j1);
k22 = sub2ind(sz, i2, j2);

% Combine linear indicies and pixel weights into one vector and sort
k      = [k11; k12; k21; k22];
w      = [w11'; w12'; w21'; w22'];
[k, I] = sort(k);
w      = w(I);

% Collapse the weights for repeated pixels
fst = [true; diff(k)~=0; true]; % marks 1st occurence of each value in k
k   = k(fst(2:end));            % unique pixel locations (linear indicies)
w   = [0; cumsum(w)];         
w   = diff(w(fst));             % collapsed weights for each pixel

% Assign the weights for all the ray sections to their pixels
d_ijm(k) = w;

% Reshape to 2D
d_ijm = reshape(d_ijm, sz);

% If needed, compute the hamming-weighted result, otherwise t_ij = 0
if obj.hamming
    t_ij = zeros((obj.Nx .^2), 1, 'double'); 

    ham_win = getWin(M, 'Hamming')';
    h11 = w11 .* ham_win;
    h12 = w12 .* ham_win;
    h21 = w21 .* ham_win;
    h22 = w22 .* ham_win;
  
    h = [h11'; h12'; h21'; h22'];
    h = h(I);

    h = [0; cumsum(h)];         
    h = diff(h(fst));              

    t_ij(k) = h;
    t_ij    = reshape(t_ij, sz);
end