%SARTEXPERIMENT Class definition for SartExperiment.
%
% DESCRIPTION:
%     SartExperiment is a class for reconstructing 2D sound speed maps from
%     time-of-flight data using straight-ray-SART (simultaneous alegbraic
%     reconstruction technique). This code is intended to be used to
%     reconstruct ultrasound tomography (UST) data acquired from 2D ring
%     array transducer hardware. 
%
% USAGE:
%     sart = SartExperiment(L, element_positions, delta_tof, temperature)
%
% INPUTS:
%     element_positions - [numeric] N x 2 matrix containing the cartesian
%                         coordinates of each transducer element. The ring
%                         array should be centred on (0, 0).
%     delta_tof         - [numeric] N x N matrix. Each (tdx, rdx) value
%                         contains the time-of-flight difference in seconds
%                         between the phantom UST data and the watershot
%                         UST data, for the ray joining the tdx and rdx
%                         transducer elements: delta_tof(tdx, rdx) =
%                         tof_phantom(tdx, rdx) - tof_water(tdx, rdx).
%                         Values of NaN should be used for missing data,
%                         not 0.
%     temperature       - [numeric] scalar water temperature during UST
%                         data acquisition [deg C]
%
% REFERENCES:
%     A. C. Kak and Malcolm Slaney, Principles of Computerized Tomographic
%     Imaging, IEEE Press, 1988 (section 7.4), available at
%     https://www.slaney.org/pct/pct-toc.html
%
% ABOUT:
%     author           - Morgan Roberts
%     date             - 9th November 2022

classdef SartExperiment < handle
    properties (GetAccess = public, SetAccess = private)
        % Do not change these values
        delta_s_ratio = 1      % length of each ray-section relative to dx
        c_vec         = [0, 0] % coordinates of array centre [m]
        
        % Properties populated in the class constructor
        detect           % structure for the detector coordinates and X/Y grid vectors
        delta_tof        % NxN input time of flight data matix
        temperature      % temperature of water during UST data acuisition [deg C]

        % Properties populated in the reconstructSart() method
        init_c_val       % scalar sound speed value used for homogeneous initial estimate [m/s]
        d                % reconstruction circle diameter [m]
        Nit              % number of iterations
        dx0              % grid step size for first iteration [m]
        upsample_factors % integer up sample factors used for each iteration (relative to dx0), must be 1 or a power of 2
        hamming          % boolean controlling whether to apply a hamming window when distributing errors along ray paths

        % Derived static properties 
        L                % physical grid length [m], automatically set to (4*obj.dx0) mm larger than reconstruction circle obj.d
        Lx               % discretised grid length [m]
        N                % Number of transducer elements in the array
        r                % reconstruction circle radius [m]
        max_Nx           % maximum number of grid points required for multi-stage reconstruction
        ray_mask         % boolean mask NxN stating which rays should be included in the reconstruction
        c_water          % sound speed of water [m/s]        

        % Derived properties that change during reconstruction
        dx               % current step size for the current iteration
        Nx               % Number of grid points (this value can change)              
        grid_x           % Nx-by-1 vector containing x and y grid points of the sound speed map 
        estimates        % sound speed estimate for each iteration [m/s]
        corrections      % weighted correction matrix for each iteration [relative slowness, s/m]
        rmses            % Room mean square error for each iteration [us]
        dxs              % grid step size for each iteration [m]
        Nxs              % grid length [pts] for each iteration
    end
    methods
        function obj = SartExperiment(element_positions, delta_tof, temperature)
            % store data
            obj.delta_tof   = delta_tof;
            obj.temperature = temperature;
            obj.c_water     = waterSoundSpeed(temperature);
            obj.ray_mask    = ~isnan(delta_tof);               

            % cartesian coordinates of detectors
            obj.detect.centroids = element_positions;
            obj.detect.x_vec     = obj.detect.centroids(:,1);
            obj.detect.y_vec     = obj.detect.centroids(:,2);   

            % compute number of transducer elements
            obj.N = size(delta_tof, 1);

        end
        
        plotSetup(obj)
        reconstructSart(obj, init_c_val, Nit, dx0, upsample_factors, options)
        plotReconResult(obj, idx, options)
        [new_c_est, weighted_correction, rmse] = performSartIteration(obj, c_est, options)
        delta_d_ijm                            = interpolateRaySection(obj, m_xvec, m_yvec, mdx)
        [m_xvec, m_yvec, M, deltaS]            = findRaySections(obj, tdx, rdx)
        saveReconResult(obj)
        
    end
end
        