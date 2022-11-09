function [m_xvec, m_yvec, M, deltaS] = findRaySections(obj, tdx, rdx)
%FINDRAYSECTIONS links a transmitter/receiver and splits the ray into sections
%
% DESCRIPTION:
%     findRaySections finds the two points where the ray between tdx and
%     rdx intersect the reconstruction circle. The line between these two
%     intersection points is divided into ray sections of length
%     obj.delta_s_ratio * obj.dx. The coordinates of the midpoint of each
%     ray section are returned, in addition to the physical length of each
%     section. If the ray doesn't intersect the circle then the output M is
%     set to zero as a flag. This function is designed to be called
%     iteratively across all tdx = 1, 2, ... N, rdx = 1, 2, ... N.
%
% USAGE:
%     [m_xvec, m_yvec, M, deltaS] = findRaySections(obj, tdx, rdx)
%
% INPUTS:
%     obj    - object, instance of the SartExperiment class
%     tdx    - [numeric] integer index selecting the current transmitter
%     rdx    - [numeric] integer index selecting the current transmitter
%
% OUTPUTS:
%     m_xvec - [numeric] M x 1 array, ray sections x-coordinates [m]
%     m_yvec - [numeric] M x 1 array, ray sections y-coordinates [m]
%     M      - [numeric] integer number of ray sections. Becomes 0 when
%              the ray does not intersect the reconstruction circle.
%     deltaS - [numeric] the physical length of each ray section [m]
%
% REFERENCES:
%     A. C. Kak and Malcolm Slaney, Principles of Computerized Tomographic
%     Imaging, IEEE Press, 1988 (section 7.4), available at
%     https://www.slaney.org/pct/pct-toc.html
%
% ABOUT:
%     author - Morgan Roberts
%     date   - 9th November 2022

% extract transmitter and receiver coordinates
r_vec = obj.detect.centroids(rdx,:)';
t_vec = obj.detect.centroids(tdx,:)';

% calculate the unit vector linking tdx/rdx
u_vec = (r_vec - t_vec) ./ ...   % unit vector giving direction
         sqrt((r_vec(1,:) - t_vec(1)).^2 + (r_vec(2,:) - t_vec(2)).^2); 

% simultaneous eqn between vector eqn of line and eqn of reconstruction circle
% coefficients for quadratic formula
a = u_vec(1,:).^2 + u_vec(2,:).^2;
b = (2 .* u_vec(1,:) .* (t_vec(1) - obj.c_vec(1))) + ...
                            (2 .* u_vec(2,:) .* (t_vec(2) - obj.c_vec(2)));
c = t_vec(1)^2 + t_vec(2)^2 + obj.c_vec(1)^2 + ...
    obj.c_vec(2)^2 - obj.r^2 - (2 * dot(t_vec, obj.c_vec));

% solutions for lambda-distance from tdx to intersection(s) with reconstruction circle
lambda1 = (-b + sqrt(b.^2 - (4 .* a .* c))) ./ (2 .* a);
lambda2 = (-b - sqrt(b.^2 - (4 .* a .* c))) ./ (2 .* a);

% cartesian coordinates of points of interception
intercept1 = t_vec + (lambda1 * u_vec);
intercept2 = t_vec + (lambda2 * u_vec);

% Physical length of intersection of ray with reconstruction circle
Li = norm(intercept1 - intercept2);

% decide number of points based on length of intersection
M = ceil(Li / (obj.dx * obj.delta_s_ratio)) + 1;

% Catch for rays with zero intersection width
if imag(lambda1) ~= 0 || M < 2
    m_xvec = [];
    m_yvec = []; 
    M = 0;
    deltaS = 0;
    return
end

% Compute physical length of each ray section
deltaS = Li / (M - 1);

% define cartesian coordinates of sections within ray
left   = min([intercept1(1), intercept2(1)]);
right  = max([intercept1(1), intercept2(1)]);
bottom = min([intercept1(2), intercept2(2)]);
top    = max([intercept1(2), intercept2(2)]);
m_xvec = linspace(left, right, M);
m_yvec = linspace(bottom, top, M);

% make sure first coordinate refers to correct circle intersection point
if m_xvec(1) ~= intercept1(1)
    m_xvec = flip(m_xvec);
end
if m_yvec(1) ~= intercept1(2)
    m_yvec = flip(m_yvec);
end

end