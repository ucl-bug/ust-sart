function plotReconResult(obj, options)
%PLOTRECONRESULT Plots the reconstruction results for the SartExperiment class
%
% DESCRIPTION:
%     plotReconResult manages plotting for the SartExperiment class. It
%     contains a subfunction plotSingleIter, which plots the sound speed
%     map, weighted correction matrix, RMSE history and step size history
%     at the specified iteration index. If the Loop option is specified,
%     plotReconResult will create an animation showing the reconstruction
%     progress.
%
% USAGE:
%     sart.plotReconResult
%     sart.plotReconResult(Iter=3, Handle=fig1, cRange=[1450,1550])
%     sart.plotReconResult(Loop=1, Handle=fig1, cRange=[1450,1550])
%
% INPUTS:
%     obj    - object, instance of the SartExperiment class
%
% OPTIONAL INPUTS:
%     Iter   - [numeric] integer specifying the iteration index at which
%              to plot the reconstruction progress. The default is the
%              final iteration. This argument is ignored if Loop is set
%              to true.
%     Handle - figure handle object to plot the result. If not specified a
%              new figure is created
%     cRange - [numeric] 2-element array in the form [c_min, c_max]
%              specifying the minimum and maximum sound speeds for the
%              colorscale. The default behaviour is autoscaling.
%     Loop   - [boolean] If set to true, the reconstruction results from
%              all iterations will be displayed on a loop.
%
% ABOUT:
%     author - Morgan Roberts
%     date   - 9th November 2022

arguments
    obj
    options.Iter   = obj.Nit;
    options.Handle = figure;
    options.cRange = [-Inf, Inf];
    options.Loop   = 0
end

if options.Loop
    for idx = 1:obj.Nit
        plotSingleIter(obj, idx, options.Handle, options.cRange);
        drawnow
        pause(0.25);
    end
else
    plotSingleIter(obj, options.Iter, options.Handle, options.cRange);
end

function plotSingleIter(obj, idx, handle, cRange)

Nx = obj.Nxs(idx);

figure(handle);
clf(handle, 'reset');
keyboard
subplot(2, 2, 1);
imagesc(1e3*obj.grid_x, 1e3*obj.grid_x, obj.estimates(1:Nx, 1:Nx, idx), cRange);
c = colorbar;
colormap(getColorMap);
title(['Iteration ', num2str(idx), ' / ', num2str(obj.Nit),]);
xlabel('x-position [mm]');
ylabel('y-position [mm]');
ylabel(c, 'Sound Speed [m/s]')
axis image

subplot(2, 2, 3);
imagesc(1e3*obj.grid_x, 1e3*obj.grid_x, obj.corrections(1:Nx, 1:Nx, idx));
c = colorbar;
colormap(getColorMap);
title('Weighted Correction');
xlabel('x-position [mm]');
ylabel('y-position [mm]');
ylabel(c, 'Relative Slowness')
axis image

subplot(2, 2, 2);
h = semilogy(1:obj.Nit, obj.rmses, 'k-');
h.Marker = '.';
h.MarkerEdgeColor = 'r';
xlabel('Iteration');
ylabel('RMSE [us]');

subplot(2, 2, 4);
h = plot(1:obj.Nit, obj.dxs*1e3, 'k-');
h.Marker = '.';
h.MarkerEdgeColor = 'r';
xlabel('Iteration');
ylabel('Grid Step Size [mm]');   

drawnow
end
end