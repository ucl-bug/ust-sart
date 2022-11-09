function plotSetup(obj)
%PLOTSETUP Plot the UST detectors, reconstruction circle and time-of-flight data
%
% DESCRIPTION:
%     plotSetup plots the UST detectors and the time of flight data matrix
%
% USAGE:
%     sart.plotSetup;
%

% % Create coordinates for the reconstruction circle
% recon_circ = makeCartCircle(obj.r, 100);
% recon_x    = [recon_circ(1,:), recon_circ(1,1)];
% recon_y    = [recon_circ(2,:), recon_circ(2, 1)];

figure;
subplot(1, 2, 1);
hold on
plot(obj.detect.x_vec, obj.detect.y_vec, 'r.');
axis image
xlabel('x-position [m]');
ylabel('y-position [m]');

subplot(1, 2, 2);
imagesc(1e6*obj.delta_tof);
c = colorbar;
axis image
xlabel('Receivers');
ylabel('Transmitters');
ylabel(c, 'Time delay [us]');
title('t_{phantom} - t_{water}');
drawnow

end