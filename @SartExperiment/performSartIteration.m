function [new_c_est, weighted_correction, rmse] = performSartIteration(obj, c_est, options)
%PERFORMSARTITERATION Manages a single SART iteration.
%
% DESCRIPTION:
%     performSartIteration takes a sound speed estimate as an input and
%     converts it to a slowness difference relative to water. This is used
%     to predict the time-of-flight at the detector positions using a
%     straight-ray model with bilinear interpolation. The difference
%     between the predicted TOF and measured TOF is then distributed to
%     pixels along the rays, according to how far the ray travelled in each
%     pixel. These correction terms are accumulated for all rays and
%     finally weighted before being added to the slowness map to give a new
%     slowness estimate. Slowness is then converted back to sound speed.
%     The RMSE time-of-flight error at the detectors is also returned.
%
% USAGE:
%     [new_c_est, weighted_correction, rmse] = performSartIteration(obj, c_est)
%     [new_c_est, weighted_correction, rmse] = performSartIteration(obj, c_est, beta=1)
%
% INPUTS:
%     obj                 - object, instance of the SartExperiment class
%     c_est               - [numeric] square matrix containing the starting
%                           sound speed estimate [m/s]
%
% OPTIONAL INPUTS:
%     beta                - [numeric] multiplication factor for the weighted 
%                           correction matrix (values greater than 1 lead
%                           to larger update steps)
%
% OUTPUTS:
%     new_c_est           - [numeric] square matrix the same size as c_est,
%                           containing the new sound speed estimate [m/s]
%     weighted_correction - [numeric] square matrix the same size as c_est,
%                           containing the weighted correction matrix
%                           [relative slowness s/m]
%     rmse                - RMSE time-of-flight error at the detectors [us]
%
% REFERENCES:
%     A. C. Kak and Malcolm Slaney, Principles of Computerized Tomographic
%     Imaging, IEEE Press, 1988 (section 7.4), available at
%     https://www.slaney.org/pct/pct-toc.html
%
% ABOUT:
%     author              - Morgan Roberts
%     date                - 9th November 2022

arguments
    obj
    c_est
    options.beta = 1;
    options.cLim = [-Inf, Inf];
end

Nx = obj.Nx;

% convert the sound speed estimate to a relative-slowness estimate
% deltaTOF_rel_water = d * relative_slowness
s_est = ((1 ./ c_est) - (1 ./ obj.ref_c));

% initialise storage
correction  = zeros(Nx, Nx);
sum_weights = zeros(Nx, Nx);
ses         = NaN * ones(Nx, Nx); % empty array for Squared detector ErrorS

[transmitters, receivers] = find(obj.ray_mask);
Nray  = nnz(obj.ray_mask); % number of rays

% Loop over rays
for idx = 1:Nray
    tdx = transmitters(idx);
    rdx = receivers(idx);  
            
    % Define the ray as a series of ray sections
    tic
    [m_xvec, m_yvec, M, deltaS] = findRaySections(obj, tdx, rdx);
    obj.ray_sect_timer = obj.ray_sect_timer + toc;

    % If ray doesn't intersect the circle, skip to next ray
    if M == 0
        continue
    end
    
    % Find contribution of each image pixel to the ray
    [d_ijm, t_ij] = calculatePixelWeights(obj, m_xvec, m_yvec);
    
    % t_ij is equal to d_ijm if hamming-weighting is turned off
    if ~obj.hamming
        t_ij = d_ijm;
    end
    obj.interpolate_timer = obj.interpolate_timer + toc;

    % calculate physical length ray spent in each pixel (deltaS is
    % the physical length of each ray section)
    a_ij = d_ijm * deltaS;

    %  calculate tof estimate for ith ray: multiply pixelwise ray distances by relative slowness
    tof_est_i = sum(a_ij .* s_est, 'all', 'omitnan'); 
    diff_i    = obj.delta_tof(tdx, rdx) - tof_est_i;

    % store the squared error in microseconds
    ses(tdx, rdx) = (diff_i * 1e6) ^ 2; 

    phys_len_i = sum(a_ij, 'all');
    % eqn 34 kak/slaney 1988 p.291
    correction = correction + (t_ij * (diff_i / phys_len_i));

    % denominator eqn 32
    sum_weights = sum_weights + d_ijm;

    % sanity check that sum equals the actual number of sections in ray
    if sum(d_ijm, 'all') - (M - 1) > 1e-3
        error('check sum matches physical length of ray')
    end

    plotting = 0;
    if plotting
        clf(gcf);
        subplot(1, 2, 1);
        imagesc(t_ij);
        title(['tx ', num2str(tdx), ' rx ', num2str(rdx)])
        axis image
        title('Ray Distance Matrix')

        subplot(1, 2, 2);
        imagesc(correction);
        title(['tx ', num2str(tdx), ' rx ', num2str(rdx)])
        axis image     
        title('Correction Matrix');
        drawnow
    end
end

% calculate root mean square error for this iteration
mse  = mean(ses, 'all', 'omitnan');
rmse = sqrt(mse); 

% calculate mask for the pixels on the reconstruction circle border
[maskX, maskY] = meshgrid(obj.grid_x, obj.grid_x);
inside = sqrt(maskX.^2 + maskY .^2) < ((obj.recon_d / 2) - obj.border_width * obj.dx0);
border = inside ~= (sum_weights > 0);

% Prevent pixels on border for being included in the update
% (Especially important for first few iterations to avoid spiking artefact)
sum_weights(border) = 0;

% calculate the weighted correction matrix
weighted_correction = correction ./ sum_weights;
% weighted_correction(isnan(weighted_correction)) = 0;
weighted_correction(sum_weights == 0) = 0;

% Calculate the new slowness estimate: apply the correction matrix from all rays
new_s_est = s_est + (options.beta * weighted_correction);

% Convert to a sound speed estimate
new_c_est = obj.ref_c ./ ((new_s_est .* obj.ref_c) + 1);

% replace Nan with zeros
new_c_est(isnan(new_c_est)) = 0;

% enforce positivity
new_c_est(new_c_est < 0) = 0;

% enforce limits
new_c_est(new_c_est < options.cLim(1)) = options.cLim(1); 
new_c_est(new_c_est > options.cLim(2)) = options.cLim(2); 

end