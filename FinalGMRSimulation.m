% FinalGMRSimulation.m

% INITIALIZE MATLAB
close all;
clc;
clear all;
set(0,'DefaultFigureWindowStyle','docked');

% UNITS
micrometers = 1;
nanometers  = 1e-3 * micrometers;

Nf = 11;
tot_ref = zeros(1,Nf);
tot_trn = tot_ref;
tot_con = tot_ref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE WAVELENGTH
lam0 = 1.55 * micrometers;
lam01 = 1.4 * micrometers;
lam02 = 1.7 * micrometers;
lam0 = linspace(lam01,lam02,Nf);

% SLAB WAVEGUIDE PARAMETERS
lamd = 1.55 * micrometers;
nslab  =  sqrt(10);
nclad1 =  1.0;
nclad2 =  sqrt(5.5);
a     = lamd/(2*nslab);
a     = 0.1583*lamd;
L    = 0.4555*lamd;             % Grating period (DO NOT CHANGE)
% d    = 0.15*lamd;                % Grating depth
d    = 0.5*a;                % Grating depth
ff   = 0.6;                     % Fill fraction  (DO NOT CHANGE)

% L    = 0.83254*lamd;             % Grating period (DO NOT CHANGE)
% d    = 0.192547096774194*lamd;             % Grating depth
% % t    = lamd/(2*nr);             % Substrate thickness
% ff   = 1-0.934096;                     % Fill fraction  (DO NOT CHANGE)

PER  = 30;                      % Number of periods
Sper = PER*L;

% GRID PARAMETERS
Sx   = (Sper + 10*L);
Sy   = 4*lamd;
NRES = 20;
NPML = [20 20 20 20];
nmax = max([nslab nclad1 nclad2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST GUESS AT RESOLUTION
dx = min([lamd lam01 lam02])/nmax/NRES;
dy = min([lamd lam01 lam02])/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
ny = ceil(a/dy);
dy = a/ny;
dx = dy;

% dx = 0.007635742857143;
% dy = 0.007727799558610;
% dx = 0.013959789473648;
% dy = 0.013937638489635;

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
ER2 = nclad1^2 * ones(Nx2,Ny2);

% CREATE SLAB
ny  = round(a/dy2);
ny1 = 1 + floor((Ny2 - ny)/2);
ny2 = ny1 + ny - 1;
nx1 = 1;
nx2 = round((PER+5)*L/dx2);
ER2(nx1:nx2,ny1:ny2) = nslab^2;
% ER2(nx2+1:end,ny1:ny2) = nslab^2;
ER2(:,ny2+1:end) = nclad2^2;

% CREATE GRATING
nd1 = ny1;
nd = round(d/dy2);
nd2 = nd1 + nd;
nd4 = ny2;
nd3 = nd4 - nd;
    nx1 = round(5*L/dx2) + round(ff*L/dx2);
ER2(1:nx1-1,nd1:nd2) = 1;
ER2(nx2+1:end,nd1:nd2) = 1;


% Chirping function
yff =@(x) x;
y   =@(x) exp(x);

for n = 1:PER
    nx1 = round((5+n)*L/dx2);
%     nx2 = nx1 + round((yff(n)/PER)*ff*lamd/dx2);
    nx2 = nx1 + round(ff*L/dx2);
    nx0 = nx1 - round((1-ff)*L/dx2);
    ER2(nx0:nx1-1,nd1:nd2-round(d*y(n/PER)/dx2)) = 1;
    ER2(nx1:nx2,nd1:nd2) = 1;
end

% SHOW ER2
figure('color','w');
imagesc(xa2./micrometers,ya2./micrometers,ER2');
xlabel('x (\mum)'); ylabel('y (\mum)');
colorbar;

figure('color','w');
lam0 = linspace(lam01,lam02,Nf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INCOPORATE PML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:Nf
clc;
disp(['Iteration ' num2str(n) ' out of ' num2str(Nf)]);
disp(['Wavelength = ' num2str(lam0(n)) ' micrometers']);
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
k0 = 2*pi/lam0(n);
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

imagesc(xa./lamd,ya./lamd,real(f)');
xlabel('x (wavelengths \lambda_0)'); ylabel('y (wavelengths \lambda_0)');
axis equal tight;
colorbar;
caxis([-1 1]*0.1)
drawnow;
fields(:,:,n) = f;
end

for a = 1:Nf
    imagesc(xa./lamd,ya./lamd,real(fields(:,:,a))');
xlabel('x (wavelengths \lambda_0)'); ylabel('y (wavelengths \lambda_0)');
axis equal tight;
colorbar;
caxis([-1 1]*0.1)
drawnow;
pause(0.1);
end