function [sx,sy] = calcpml2d(NGRID,NPML)
% CALCPML2D Calculate the PML parameters on a 2D grid
%
% [sx,sy] = calcpml2d(NGRID,NPML);
%
% This MATLAB function calculates the PML parameters sx and sy
% to absorb outgoing waves on a 2D grid.
%
% Input Arguments
% =================
% NGRID Array containing the number of points in the grid
% = [ Nx Ny ]
% NPML Array containing the size of the PML at each boundary
% = [ Nxlo Nxhi Nylo Nyhi ]
%
% Output Arguments
% =================
% sx,sy 2D arrays containing the PML parameters on a 2D grid

% CONSTANTS
n0 = 376.73;

% INITIALIZE VARIABLES
Nx      = NGRID(1);
Ny      = NGRID(2);
NxLo    = NPML(1);
NxHi    = NPML(2);
NyLo    = NPML(3);
NyHi    = NPML(4);
amax    = 3;
p       = 3;
sigmax  = 1;

% INITIALIZE sx AND sy
sx = ones(Nx,Ny);
sy = sx;

% ADD NxLo PML
for nx = 1 : NxLo
    ax = 1 + amax *(nx/NxLo)^p;
    sigx = sigmax *(sin(0.5*pi*nx/NxLo))^2;
    sx(NxLo-nx+1,:) = ax * (1 + (n0*sigx) * 1i);
end

% ADD NxHi PML
for nx = 1 : NxHi
    ax = 1 + amax *(nx/NxHi)^p;
    sigx = sigmax *(sin(0.5*pi*nx/NxHi))^2;
    sx(Nx-NxHi+nx,:) = ax * (1 + (n0 * sigx) * 1i);
end

% ADD NyLo PML
for ny = 1 : NyLo
    ay = 1 + amax * (ny/NyLo)^p;
    sigy = sigmax * (sin(0.5*pi*ny/NyLo))^2;
    sy(:,NyLo-ny+1) = ay * (1 + (n0 * sigy) * 1i);
end

% ADD NyHi PML
for ny = 1 : NyHi
    ay = 1 + amax * (ny/NyHi)^p;
    sigy = sigmax * (sin(0.5*pi*ny/NyHi))^2;
    sy(:,Ny-NyHi+ny) = ay * (1 + (n0 * sigy) * 1i);
end

end