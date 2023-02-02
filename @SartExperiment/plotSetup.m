function plotSetup(obj, options)
%PLOTSETUP Plot the UST detectors, reconstruction circle and time-of-flight data
%
% DESCRIPTION:
%     plotSetup plots the UST detectors and the time of flight data matrix
%
% USAGE:
%     sart.plotSetup;
%     sart.plotSetup(recon_d=recon_d);
%
% INPUTS:
%     obj           - object, instance of the SartExperiment class
%
% OPTIONAL INPUTS:
%     recon_d       - [numeric] scalar diameter of the reconstruction circle
%                     to plot [m]. If not supplied, uses the default value
%                     calculated in the SartExperiment class constructor.
%
% ABOUT:
%     author        - Morgan Roberts
%     date          - 2nd February 2023

arguments
    obj
    options.recon_d = obj.default_recon_d;
end

% Create coordinates for the reconstruction circle
recon_circ = makeCartCircle(options.recon_d/2, 100);
recon_x    = [recon_circ(1,:), recon_circ(1,1)];
recon_y    = [recon_circ(2,:), recon_circ(2,1)];

figure;
subplot(1, 2, 1);
hold on
plot(obj.detect.x_vec, obj.detect.y_vec, 'r.');
plot(recon_x, recon_y, 'k-');
axis image
xlabel('x-position [m]');
ylabel('y-position [m]');
legend({'Detectors', 'Reconstruction Circle'}, location='northoutside');

subplot(1, 2, 2);
imagesc(1e6*obj.delta_tof);
c = colorbar(location='northoutside');
axis image
xlabel('Receivers');
ylabel('Transmitters');
ylabel(c, 'Time delay [us]');
title('t_{phantom} - t_{water}');
drawnow

end