% time_of_arrival_example.m
%
% DESCRIPTION:
%     This script is a tutorial with three examples that demonstrate the
%     use of the pickFirstMotionAIC function for simulated and experimental
%     data.
%
% ABOUT
%     author    - Morgan Roberts
%     date      - 22/11/2022

close all
clearvars

load('time_of_arrival_example_data.mat', 'sim', 'expA', 'expB');

% -------------------------------------------------------------------------
% Example 1: Simulated Data

% It fails due to the leading zeros
pickFirstMotionAIC(sim, Plot=true);

% Try again with some added noise - it succeeds
i_fm = pickFirstMotionAIC(sim, Plot=true, NoiseLevel=0.5e-2);
disp(['Example 1, first-motion-sample: ', num2str(i_fm)]);

% -------------------------------------------------------------------------
% Example 2: Experimental Data A

% It has an early-picking-error due to artefact before the first arrival
pickFirstMotionAIC(expA, Plot=true);

% Try again with some added noise - it succeeds
i_fm = pickFirstMotionAIC(expA, Plot=true, NoiseLevel=1.5e-2);
disp(['Example 2, first-motion-sample: ', num2str(i_fm)]);

% -------------------------------------------------------------------------
% Example 3: Experimental Data B

% It has a late-picking-error
pickFirstMotionAIC(expB, Plot=true, NoiseLevel=1.5e-2);

% Try again, but ignore the ring-down after the main peak - it succeeds
i_fm = pickFirstMotionAIC(expB, Plot=true, NoiseLevel=1.5e-2, IgnoreRingDown=true);
disp(['Example 3, first-motion-sample: ', num2str(i_fm)]);

% You could instead pre-pad the signal with zeros first
expB = [zeros(1, 50), expB];
pickFirstMotionAIC(expB, Plot=true, NoiseLevel=1.5e-2, IgnoreRingDown=false);