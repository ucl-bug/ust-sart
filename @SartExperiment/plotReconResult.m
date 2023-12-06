function plotReconResult(obj, customAxes, options)
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
    customAxes
    options.Iter   = obj.Nit;
    options.cRange = [-Inf, Inf];
    options.Loop   = 0;
end

if options.Loop
    for idx = 1:obj.Nit
        plotSingleIter(obj, idx, options.cRange, customAxes);
        drawnow
        pause(0.25);
    end
else
    plotSingleIter(obj, options.Iter, options.cRange, customAxes);
end

function plotSingleIter(obj, idx, cRange, customAxes)

Nx = obj.Nxs(idx);

% figure(handle);
% clf(handle, 'reset');

% subplot(, 2, 2, 1);
% data = obj.estimates(1:Nx, 1:Nx, idx);
% maxDev = max(abs([max(data(:)) - obj.ref_c, min(data(:)) - obj.ref_c]));
% cRange = [obj.ref_c - maxDev, obj.ref_c + maxDev];

imagesc(customAxes(1), obj.estimates(1:Nx, 1:Nx, idx), cRange);
c = colorbar(customAxes(1), 'southoutside');
colormap(customAxes(1), getColorMap);
colormap(customAxes(1), 'gray');
title(customAxes(1), ['Iteration ', num2str(idx), ' / ', num2str(obj.Nit),]);
xlim(customAxes(1), [1, Nx]);
ylim(customAxes(1), [1, Nx]);
customAxes(1).XAxis.Visible = 'off';
customAxes(1).YAxis.Visible = 'off';
ylabel(c, '[m/s]');
axis(customAxes(1), 'square');

% subplot(2, 2, 3);
imagesc(customAxes(3), obj.corrections(1:Nx, 1:Nx, idx));
colormap(customAxes(3), getColorMap);
colormap(customAxes(3), 'gray');
xlim(customAxes(3), [1, Nx]);
ylim(customAxes(3), [1, Nx]);
title(customAxes(3), 'Weighted Correction');
customAxes(3).XAxis.Visible = 'off';
customAxes(3).YAxis.Visible = 'off';
axis(customAxes(3), 'square');

% subplot(, 2, 2, 2);
h = semilogy(customAxes(2), 1:obj.Nit, obj.rmses, 'k-');
h.Marker = '.';
h.MarkerEdgeColor = 'r';
xlabel(customAxes(2), 'Iteration');
ylabel(customAxes(2), 'RMSE [\mus]');
customAxes(2).XAxis.Visible = 'off';
title(customAxes(2), '');

% subplot(, 2, 2, 4);
h = plot(customAxes(4), 1:obj.Nit, obj.dxs*1e3, 'k-');
h.Marker = '.';
h.MarkerEdgeColor = 'r';
xlabel(customAxes(4), 'Iteration');
ylabel(customAxes(4), 'dx [mm]');   
title(customAxes(4), '');

drawnow
end
end