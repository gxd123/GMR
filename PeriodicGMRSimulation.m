% FinalGMRSimulation.m

% INITIALIZE MATLAB
close all;
clc;
clear all;
set(0,'DefaultFigureWindowStyle','docked');

% UNITS
micrometers = 1;
nanometers  = 1e-3 * micrometers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE WAVELENGTH
lam0 = 1.55 * micrometers;

% SLAB WAVEGUIDE PARAMETERS
lamd = 1.5605 * micrometers;
nslab = sqrt(2.0);
nclad = 1.0;
a     = lamd/(2*nslab);
L    = 0.8621*lamd;             % Grating period (DO NOT CHANGE)
d    = 0.15*lamd;                % Grating depth
ff   = 0.5;                     % Fill fraction  (DO NOT CHANGE)
PER  = 200;                      % Number of periods

% GRID PARAMETERS
Sx   = 210*lamd;
Sy   = 10*lamd;
NRES = 20;
NPML = [20 20 20 20];
nmax = max([nslab nclad]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST GUESS AT RESOLUTION
dx = lamd/nmax/NRES;
dy = lamd/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
ny = ceil(a/dy);
dy = a/ny;
dx = dy;
% dx = 0.007635742857143;
% dy = 0.007727799558610;

% GRID SIZE
Nx = NPML(1) + ceil(Sx/dx) + NPML(2);
Sx = Nx*dx;
Ny = NPML(3) + ceil(Sy/dy) + NPML(4);
Sy = Ny*dy;
% Ny = Nx;

% 2X GRID
Nx2 = 2*Nx;     dx2 = dx/2;
Ny2 = 2*Ny;     dy2 = dy/2;

% GRID AXES
xa = [0:Nx-1]*dx;
ya = [0:Ny-1]*dy;
xa2 = [0:Nx2-1]*dx2;
ya2 = [0:Ny2-1]*dy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATERIALS
UR2 = ones(Nx2,Ny2);
ER2 = nclad^2 * ones(Nx2,Ny2);

% CREATE SLAB
ny  = round(a/dy2);
ny1 = 1 + floor((Ny2 - ny)/2);
ny2 = ny1 + ny - 1;
nx1 = 1;
nx2 = round((PER+5)*lamd/dx2);
ER2(nx1:nx2,ny1:ny2) = nslab^2;

% CREATE GRATING
nd1 = ny1;
nd = round(d/dy2);
nd2 = nd1 + nd;
nd4 = ny2;
nd3 = nd4 - nd;

% Chirping function
yff =@(x) 0;
y   =@(x) 0;

for n = 1:PER
    nx1 = round((5+n)*lamd/dx2);
%     nx2 = nx1 + round((yff(PER-n)/PER)*ff*lamd/dx2);
    nx2 = nx1 + round(ff*lamd/dx2);
    ER2(nx1:nx2,nd1:nd2-round(d*y((PER-n)/PER)/dx2)) = 1;
%     ER2(nx1:nx2,nd3:nd4) = 1;
end

% SHOW ER2
figure('color','w');
imagesc(xa2./lamd,ya2./lamd,ER2');
xlabel('x (wavelengths \lambda_0)'); ylabel('y (wavelengths \lambda_0)');
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INCOPORATE PML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE PML PARAMETERS
[sx,sy] = calcpml2d([Nx2 Ny2],2*NPML);

% INCORPORATE PML
ERxx = ER2 ./ sx .* sy;
ERyy = ER2 .* sx ./ sy;
ERzz = ER2 .* sx .* sy;

URxx = UR2 ./ sx .* sy;
URyy = UR2 .* sx ./ sy;
URzz = UR2 .* sx .* sy;

% PARSE TO 1X GRID
ERxx = ERxx(2:2:Nx2,1:2:Ny2);
ERyy = ERyy(1:2:Nx2,2:2:Ny2);
ERzz = ERzz(1:2:Nx2,1:2:Ny2);

URxx = URxx(1:2:Nx2,2:2:Ny2);
URyy = URyy(2:2:Nx2,1:2:Ny2);
URzz = URzz(2:2:Nx2,2:2:Ny2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT 1D CROSS SECTION
nxs  = NPML(1) + 2;
ny1  = NPML(3) + 1;
ny2  = Ny - NPML(4);
erzz = ERzz(nxs,ny1:ny2);
urxx = URxx(nxs,ny1:ny2);
uryy = URyy(nxs,ny1:ny2);

% DIAGONALIZE MATERIALS
erzz = diag(sparse(erzz));
urxx = diag(sparse(urxx));
uryy = diag(sparse(uryy));

% BUILD DERIVATE MATRICES
k0 = 2*pi/lam0;
ny = ny2 - ny1 + 1;
[DEX,DEY,DHX,DHY] = yeeder([1 ny],k0*[dx dy],[0 0 0 0]);

% BUILD EIGEN-VALUE PROBLEM
A = - (DHY/urxx*DEY + erzz);
B = inv(uryy);

% SOLVE EIGEN-VALUE PROBLEM
[V,D] = eig(full(A),full(B));
D = diag(D);

% IDENTIFY FUNDAMENTAL MODE
Ez0 = V(:,1);
neff = sqrt(-D(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDFD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BUILD DERIVATE MATRICES
[DEX,DEY,DHX,DHY] = yeeder([Nx Ny],k0*[dx dy],[0 0 0 0]);

% DIAGONALIZE MATERIALS
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% BUILD WAVE MATRIX
A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;

% SOURCE FIELD
fsrc = zeros(Nx,Ny);
for nx = 1 : Nx
    fsrc(nx,ny1:ny2) = Ez0*exp(1i*k0*neff*nx*dx); 
end

% Q
Q = zeros(Nx,Ny);
Q(1:nxs,:) = 1;
Q = diag(sparse(Q(:)));

% SOURCE VECTOR
b = (Q*A - A*Q)*fsrc(:);

% CALCULATE FIELD
f = A\b;
f = reshape(f,Nx,Ny);

figure('color','w');
imagesc(xa./lamd,ya./lamd,real(f)');
xlabel('x (wavelengths \lambda_0)'); ylabel('y (wavelengths \lambda_0)');
axis equal tight;
colorbar;