function [gridVec, dx, soundSpeed] = getFinalEstimate(obj)

% Find the last iteration index to successfully complete
idx = obj.completedUpto;

% Extract the grid size and length
dx = obj.dxs(idx);
Nx = obj.Nxs(idx);

% Work out what the grid vector was
gridVec = (0:(Nx-1)) * dx;
gridVec = gridVec - gridVec(ceil(Nx / 2));

% Extract the sound speed
soundSpeed = obj.estimates(1:Nx, 1:Nx, idx);

end


