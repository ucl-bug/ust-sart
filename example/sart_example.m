close all
clearvars

% Make sure ust-sart is on the path
addpath(genpath('..\')); 
load('example_element_positions');
load('example_delta_tof');

% remove lower triangle TOF data (it's reciprocal) to save computation time
delta_tof(logical(tril(ones(size(delta_tof)), -1))) = NaN;

% Initialise sart object
L           = 0.23; % grid length [m]
temperature = 20;   % water temperature [degC]
sart        = SartExperiment(L, element_positions, delta_tof, temperature);

% Plot setup and data
sart.plotSetup;

% Perform reconstruction
ups        = [1, 1];       % upsampling factors
Nit        = length(ups);  % number of iterations
dx0        = 4e-3;         % step size for iteration 1 [m]
init_c_val = sart.c_water; % sound speed value for homogeneous initial estimate [m/s]
hamming    = 0;            % boolean controlling whether hamming window is used
recon_d    = 0.135;        % diameter of reconstruction circle [m]
sart.reconstructSart(init_c_val, recon_d, Nit, dx0, ups, hamming=hamming);

% plot final estimate
figure; imagesc(squeeze(sart.estimates(:,:,end))); axis image;