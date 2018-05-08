%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Resonance Grating - FDFD Simulation
%       This MATLAB files simulates a single period of a rectangular
%       grating using the Finite-Difference Frequency Domain method.
%       The wavelength sweep goes from 2 um to 1 um.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'DefaultFigureWindowStyle','docked');

% close all;
clc;
clear all;
close all;

% UNITS
micrometers = 1;
millimeters = 1e3 * micrometers;
meters = 1e3 * millimeters;
degrees = pi/180;
seconds = 1;
hertz = 1/seconds;
gigahertz = 1e9 * hertz;

% CONSTANTS
c0 = 299792458 * meters/seconds;

% POINTS FOR SWEEP 
Nf = 101;
tot_ref = zeros(1,Nf);
tot_trn = tot_ref;
tot_con = tot_ref;

% FIGURE SETTINGS 
fig = 1;        % 0 for no figures, 1 for figure animation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE PARAMETERS
lam0 = 1.55 * micrometers;
f0 = c0/lam0;
lam01 = 1.0 * micrometers;
lam02 = 2.1 * micrometers;
lam0f = linspace(lam01,lam02,Nf);

% GRATING PARAMETERS
lamd   = 1.55 * micrometers;        % Design wavelength
fd     = c0/lamd;                   % Design frequency
nclad1 = 1.0;                       % Reflection Region
nclad2 = 1.5;                       % Transmission Region
nslab  = 2.0;                       % Waveguide
ngrat  = 1.7;                       % Grating effective refractive index
L      = 0.5000*lamd;               % Grating period
ff     = 0.5;                       % Grating fill fraction
d      = 0.1500*lamd;               % Grating depth
t      = lamd/(4*nslab);            % Substrate thickness

% EXTERNAL MATERIALS
ur1 = 1.0;                    % Reflection region permeability
er1 = nclad1^2;               % Reflection region permittivity
ur2 = 1.0;                    % Transmission region permeability
er2 = nclad2^2;               % Transmission region permittivity

% GRID PARAMETERS
NRES = 60;                    % Grid resolution
BUFZ = 2*lam02 * [1 1];       % Spacer region above and below grating
DEV.NPML = [20 20];           % Size of PML at top and bottom of grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider wavelengths
lam_min = min([lam01 lam02])/max([nclad1 nclad2 nslab]);
dlam = lam_min/NRES; 

% Consider mechanical parameters
dmin = ff*L;   % Fill fraction is
dd = dmin/2;   % Delta for distance

% Choose the highest resolution
dx = dlam;
dy = dx;

% dx =    0.055172006602080
% dy = dx;

% Snap grid to critical dimension (in this case L and d+t)
Nx = 2*ceil(L/dx/2) + 1;      % First guess at grid spacing (odd for periodic)
Ny = ceil(d/dy);

% Calculate new grid resolutions
dx = L/Nx;
dy = d/Ny;

% Incorporate PML and spacer regions
Ny = Ny + DEV.NPML(1) + DEV.NPML(2) + ceil(BUFZ(1)/dy) + ceil(BUFZ(2)/dy);

% Calculate 2x grid
Nx2 = 2*Nx;
Ny2 = 2*Ny;
dx2 = dx/2;
dy2 = dy/2;
xa2 = [0:Nx2-1]*dx2;
ya2 = [0:Ny2-1]*dy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON 2x GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UR2 = ur1*ones(Nx2,Ny2);
ER2 = er1*ones(Nx2,Ny2);

% Calculate start and stop indices for filling in the grid
nt1 = 2*DEV.NPML(1) + ceil(BUFZ(1)/dy2);
nt2 = nt1 + round(t/dy2);
nt3 = nt2 + 1;
nt4 = nt3 + round(d/dy2);
nd2 = nt2;
nd1 = nt2 - round(d/dy2);
nx1 = round(Nx2*ff/2);
nx2 = round(Nx2 - Nx2*ff/2);

% Fill in the permeability regions
UR2(:,:) = 1.0;

% Fill in the permittivity regions
ER2(:,nt1:nt2) = nslab^2;
ER2(:,nt3:nt4) = ngrat^2;
ER2(:,1:nt1-1) = er2;
% ER2(nx1:nx2,nd1:nd2) = 1;

ER2 = fliplr(ER2);
figure('color','w');
if fig
    subplot(2,2,1:2);
end
imagesc(xa2,ya2,ER2);
title('\epsilon_r');
xlabel('x (\mum)'); ylabel('y (\mum)');
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TMM 3D SWEEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf = linspace(0.01*lamd,lamd,Nf);
df = linspace(0.01*lamd,lamd,Nf);
nf = linspace(1,nslab,Nf);

% Device parameters
DEV.er1 = er1;
DEV.er2 = er2;
DEV.ur1 = ur1;
DEV.ur2 = ur2;
DEV.ER  = [ngrat^2 nslab^2];
DEV.UR  = [1.0 1.0];
DEV.L   = [d t];
DEV.NP  = 1;

% Source parameters
SRC.lam0  = 1.55 * micrometers;    % Free space wavelength
SRC.theta = 0 * degrees;
SRC.phi   = 0 * degrees;
SRC.ate   = 1;
SRC.atm   = 0;

for n = 1:Nf                % Sweep grating depth
    tic;
    for m = 1:Nf            % Sweep slab thickness
        for k = 1:Nf        % Sweep grating effective refractive index
            DEV.L  = [df(n) tf(m)];
            DEV.ER = [nf(k)^2 nslab^2];
            DAT = tmm1d(DEV,SRC);
            tot_ref(n,m,k) = DAT.REF;
%             tot_trn(n,m,k) = DAT.TRN;
%             tot_con(n,m,k) = tot_ref(n) + tot_trn(n);
            

%             subplot(2,2,3:4);
%             plot(df(1:n)./micrometers,100.*tot_ref(1:n),'r','linewidth',2);
%             hold on;
%             plot(df(1:n)./micrometers,100.*tot_trn(1:n),'b','linewidth',2);
%             plot(df(1:n)./micrometers,100.*tot_con(1:n),'--k','linewidth',2);
%             hold off;
%             title(['Grating Thickness Sweep']);
%             xlabel('Grating Depth(\mum)'); ylabel('Power (%)');
%             legend('Reflectance','Transmittance','Conservation');
%             xlim([df(1) df(end)]./micrometers); ylim([0 102]);
%             drawnow;
        end
    end
    clc;
            time = toc;
            mint = floor(time*(Nf-n)/60);
            sec = round((time*(Nf-n)/60 - mint)*60);
            disp(['Estimated Time: ' num2str(mint) ' minutes '...
                    num2str(sec) ' seconds']);
            disp([num2str(n) ' out of ' num2str(Nf) ' Iterations']);
end

figure('color','w');
for a = 1:Nf
surf(df./micrometers,tf./micrometers,tot_ref(:,:,a))
title(['Reflectance Effective Refractive Index ' num2str(nf(a))]);
xlabel('Grating Depth (\mum)');
ylabel('Slab Thickness (\mum)');
zlabel('Reflectance');
drawnow;
pause(1);
end
