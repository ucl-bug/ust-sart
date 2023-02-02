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

% Visualise default reconstruction circle diameter, calculated for this geometry
recon_d = sart.default_recon_d;
sart.plotSetup(recon_d=recon_d);

% Perform reconstruction
ups          = [1 * ones(1, 45), ...
                2 * ones(1, 30),...
                4 * ones(1, 5)]; % upsampling factors for each iteration
Nit          = length(ups);  % number of iterations
dx0          = 4e-3;         % step size for iteration 1 [m]
init_c_val   = sart.c_water; % sound speed value for homogeneous initial estimate [m/s]
hamming      = 1;            % boolean controlling whether hamming window is used
border_width = 2;            % width of non-updating region at edge of reconstruction circle (integer multiples of dx0)
sart.reconstructSart(init_c_val, Nit, dx0, ups, ...
    recon_d=0.135, border_width=border_width, hamming=hamming);

% plot final estimate and save
sart.plotReconResult(cRange=[1425, 1580]);
sart.saveReconResult;