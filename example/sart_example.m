close all
clearvars

% Make sure ust-sart is on the path
addpath(genpath('..\')); 
load('example_element_positions');
load('example_delta_tof');

% remove lower triangle TOF data (it's reciprocal) to save computation time
delta_tof(logical(tril(ones(size(delta_tof)), -1))) = NaN;

% Initialise sart object
sart = SartExperiment(element_positions, delta_tof);

% Visualise default reconstruction circle diameter, calculated for this geometry
recon_d = sart.default_recon_d;
sart.plotSetup(recon_d=recon_d);

% Perform reconstruction
ups          = [1 * ones(1, 45), ...
                2 * ones(1, 30),...
                3 * ones(1, 5)]; % upsampling factors for each iteration
dx0          = 4e-3;         % step size for iteration 1 [m]
ref_c        = 1480;         % reference sound speed value for background [m/s]
hamming      = true;         % boolean controlling whether hamming window is used
border_width = 1;            % width of non-updating region at edge of reconstruction circle (integer multiples of dx0)
sart.reconstructSart(ref_c, dx0, ups, ...
    recon_d=0.135, border_width=border_width, hamming=hamming);

% plot final estimate and save
sart.plotReconResult(cRange=[1425, 1580]);
sart.saveReconResult;