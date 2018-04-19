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

% PML CHARACTERISTICS
n0 = 376.73031346177;   % Free space impedance
amax = 3;               % Amplitude term
p = 3;                  % Exponential term
omax = 1;               % Conductivity term

% INITIALIZE ALL PML MATRICES TO ONES
sx = ones(NGRID(1),NGRID(2));
sy = ones(NGRID(1),NGRID(2));

% FILL IN X-AXIS PML REGIONS
% ADD XLO PML
for nx = 1 : NPML(1)
    sx(NPML(1)-nx+1,:) = (1+amax*(nx/NPML(1))^p) *...
        (1+1i*n0*omax*sin(pi/2 * nx/NPML(1))^2);
end

% ADD XHI PML
for nx = 1 : NPML(2)
    sx(NGRID(1)-NPML(2)+nx,:) = (1+amax*(nx/NPML(2))^p) *...
        (1+1i*n0*omax*sin(pi/2 * nx/NPML(2))^2);
end

% ADD YLO PML
for ny = 1 : NPML(3)
    sy(:,NPML(3)-ny+1) = (1+amax*(ny/NPML(3))^p) *...
        (1+1i*n0*omax*sin(pi/2 * ny/NPML(3))^2);
end

% ADD YHI PML
for ny = 1 : NPML(4)
    sy(:,NGRID(2)-NPML(4)+ny) = (1+amax*(ny/NPML(4))^p) *...
        (1+1i*n0*omax*sin(pi/2 * ny/NPML(4))^2);
end


