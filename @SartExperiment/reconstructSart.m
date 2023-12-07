function reconstructSart(obj, ref_c, dx0, upsample_factors, options)
%RECONSTRUCTSART runs an interative SART reconstruction for sound speed inversion of time-of-flight data.
%
% DESCRIPTION:
%     reconstructSart manages the iterations for a SART experiment. A
%     homogeneous initial sound speed estimate is iteratively updated using
%     Kaczmarz's method of projections to find the sound speed map that
%     minimises the error between simulated and measured time-of-flight
%     data at the detector positions. reconstructSart stores the RMSE
%     error, the correction matrix, and the sound speed estimate, for every
%     iteration. Each iteration step is also plotted.
%
% USAGE:
%     reconstructSart(obj, init_c_val, Nit, dx0, upsample_factors)
%     reconstructSart(init_c_val, Nit, dx0, upsample_factors, Hamming=1)
%     
% INPUTS:
%     obj              - object, instance of the SartExperiment class
%     ref_c            - [numeric] Scalar reference (background) sound speed value [m/s].
%     Nit              - [numeric] Integer number of iterations to perform
%     dx0              - [numeric] Grid step size for the first iteration [m]
%     upsample_factors - [numeric] 1 x Nit array containing the up-sample
%                        factor to use at each iteration (relative to dx0).
%                        Must monotonically increase. 
%                        Each element must be 1, or a power of 2.
%
% OPTIONAL INPUTS:
%     hamming          - [boolean] Set to true to use a Hamming window when
%                        distributing errors along each ray
%                        (t_ij in eqn 35 section 7.4 in Kak&Slaney1988)
%     recon_d          - [numeric] scalar diameter of the reconstruction
%                        circle [m]. If not supplied, uses the default
%                        value calculated in the SartExperiment class
%                        constructor.
%     border_width     - [numeric] Integer number of dx0 grid points that
%                        are prevented from updating at the edge of the
%                        reconstruction circle to avoid spiking errors.
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
    ref_c
    dx0
    upsample_factors

    options.hamming = 0;
    options.recon_d = obj.default_recon_d;
    options.border_width = 1;
    options.customAxes = [];
    options.cLim = [-Inf, Inf];
end

border_proportion = (2 * (options.border_width * dx0)) / options.recon_d;
if border_proportion > 0.15
    warning(['The border is currently ', num2str(border_proportion*1e2), ...
        ' % of the reconstruction circle. Pixels in the border are not',...
        ' updated. Consider decreasing border_width, or increasing dx0.']);
    fprintf('\n');
end

if options.recon_d > obj.default_recon_d
    warning(['Reconstruction circle diameter ', num2str(1e3*options.recon_d), ...
        ' mm is large. The circle should cover the expected diameter of', ...
        ' the reconstructed object only. Reduce this value to less than', ...
        num2str(1e3*obj.default_recon_d), ...
        ' mm to increase reconstruction speed']);
end

% -----------------------
% SETUP
% -----------------------

% Store values
obj.recon_d          = options.recon_d;
obj.r                = obj.recon_d / 2;
obj.hamming          = options.hamming;
obj.dx0              = dx0;
obj.upsample_factors = upsample_factors;
obj.ref_c            = ref_c;
obj.Nit              = length(upsample_factors);
obj.border_width     = options.border_width;

% set the physical grid length [m]
pad   = 2; % grid points to pad outside reconstruction circle
obj.L = obj.recon_d + (obj.dx0 * pad * 2);

% Round the grid length to fit an odd integer number of grid points
obj.Nx = 2 * ceil(obj.L / (2 * obj.dx0)) + 1;

% calculate the actual discretised grid length [m] - constant throughout
obj.Lx = (obj.Nx - 1) * obj.dx0;

% recalculate the starting grid size (this changes due to upsampling)
obj.Nx = 2 * ceil(obj.Lx / (2 * obj.dx0)) + 1;

% calculate maximum number of grid points required across all iterations
max_up     = max(obj.upsample_factors);
min_dx     = obj.dx0 / max_up;
obj.max_Nx = 2 * ceil(obj.Lx / (2 * min_dx)) + 1;

% initialise memory arrays
obj.rmses       = NaN * ones(1, obj.Nit);
obj.estimates   = zeros(obj.max_Nx, obj.max_Nx, obj.Nit);
obj.corrections = zeros(obj.max_Nx, obj.max_Nx, obj.Nit);
obj.dxs         = NaN * ones(1, obj.Nit);
obj.Nxs         = NaN * ones(1, obj.Nit);

% create the grid vector for the initial estimate sound speed map
obj.grid_x = (0:(obj.Nx-1)) * obj.dx0;
obj.grid_x = obj.grid_x - obj.grid_x(ceil(obj.Nx / 2));

% -----------------------
% RUN RECONSTRUCTION
% -----------------------
disp('----------------------------------------------------------');
if isempty(options.customAxes)
    figure('Position', [576 405 761 411]);
    ax1 = subplot(3, 2, [1, 3, 5]);
    ax2 = subplot(3, 2, 4);
    ax3 = subplot(3, 2, 2);
    ax4 = subplot(3, 2, 6);
    options.customAxes = [ax1, ax2, ax3, ax4];
end

for idx = 1:obj.Nit

    % set the grid step size for this iteration and store
    obj.dx       = obj.dx0 / obj.upsample_factors(idx);
    obj.dxs(idx) = obj.dx;

    % extract the starting sound-speed estimate for this iteration
    if idx == 1
        c_est = obj.ref_c * ones(obj.Nx);
    else
        c_est = obj.estimates(1:obj.Nx, 1:obj.Nx, idx-1);
    end

    % update the number of grid points required and upsample the grid vector
    obj.Nx    = 2 * ceil(obj.Lx / (2 * obj.dx)) + 1;
    grid_x_up = (0:(obj.Nx-1)) * obj.dx;
    grid_x_up = grid_x_up - grid_x_up(ceil(obj.Nx / 2));
    obj.Nxs(idx) = obj.Nx;

    t_iter = tic;
    disp(['Iter ', num2str(idx), ' / ', num2str(obj.Nit), ' dx = ', num2str(1e3*obj.dx), ' mm, Nx = ', num2str(obj.Nx), ' pts']);
    disp('Reconstructing...');
    
    % interpolate the sound speed estimate to match the new grid vector
    if ~isequal(size(c_est), [obj.Nx, obj.Nx])
        c_est = interp2(obj.grid_x, obj.grid_x, c_est, grid_x_up, grid_x_up', 'cubic');
    end

    % activate the new upsampled grid vector
    obj.grid_x = grid_x_up;
    
    % calculate new sound speed estimate
    [new_c_est, weighted_correction, rmse] = performSartIteration(obj, c_est, beta=1, cLim=options.cLim);

    % store estimates, weighted correction and RMSE
    obj.rmses(idx) = rmse;
    obj.estimates(1:obj.Nx, 1:obj.Nx, idx) = new_c_est;
    obj.corrections(1:obj.Nx, 1:obj.Nx, idx) = weighted_correction;

    % Update the plot
    plotReconResult(obj, options.customAxes, Iter=idx)

    % Display iteration time stat
    obj.total_timer = obj.total_timer + toc(t_iter);
    disp(['Completed in ', num2str(toc(t_iter))]);   

    % indicate that this iteration finished without errors or interuption
    obj.completedUpto = idx;

end

disp('----------------------------------------------------------');
disp(['Time for breaking rays into sections: ', num2str(obj.ray_sect_timer), ' s']);
disp(['Time for calculating pixel weights:   ', num2str(obj.interpolate_timer), ' s']);
disp(['Total time:                           ', num2str(obj.total_timer), ' s']);
disp('----------------------------------------------------------');

end