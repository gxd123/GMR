%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMR_Resonance.m
% This script simulates a 1D GMR Filter to find its resonant wavelength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESTORE STATE
clear all;
close all;
clc;

% OPEN FIGURE WINDOW
fig  = figure('Color','w');

% UNITS
nanometers  = 1;
micrometers = 1e3 * nanometers;
millimeters = 1e3 * micrometers;
meters      = 1e3 * millimeters;
degrees     = pi/180;
seconds     = 1;
hertz       = 1/seconds;
gigahertz   = 1e9 * hertz;

% CONSTANTS
c0 = 299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
NLAM        = 501;
lam1        = 400 * nanometers;
lam2        = 700 * nanometers;
lam0        = linspace(lam1,lam2,NLAM);
SRC.theta   = 0 * degrees;
SRC.MODE    = 'E';

% GRATING PARAMETERS
T   = 134 * nanometers; % Depth of deposition
LAM = 314 * nanometers; % Grating period
x1  = 0.5 * LAM;        % Fill fraction of material 1
x2  = 0.5 * LAM;        % Fill fraction of material 2
nL  = 2.0;              % Low n
eL  = nL^2;             % Permittivity of low index material
uL  = 1;                % Permeabillity of low index material
nH  = 2.1;              % High n
eH  = nH^2;
uH  = 1;

% EXTERNAL MATERIALS
n1  = 1.0;
ur1 = 1.0; %permeability in the reflection region
er1 = n1^2; %permittivity in the reflection region
n2  = 1.52;
ur2 = 1.0; %permeability in the transmission region
er2 = n2^2; %permittivity in the transmission region

% GRID PARAMETERS
NRES = 100;
BUFZ = 4*max(lam0) * [1 1];
DEV.NPML = [0 0 20 20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE REFRACTIVE INDICES
nmax    = max([n1 n2 nL nH]);

% RESOLVE FOR MINIMUM WAVELENGTH
lam_min = min(lam0);

% COMPUTE INITIAL GRID RESOLUTION
dx      = lam_min/nmax/NRES;
dy      = lam_min/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSION
Nx      = 2*ceil(LAM/dx/2) + 1; % Ensure grid is always odd
dx      = LAM/Nx;
Ny      = ceil(T/dy); 
dy      = T/Ny; % Recalculate resolution

% COMPUTE GRID SIZE
Sx      = LAM;
Sy      = T + sum(BUFZ);
Ny      = ceil(Sy/dy) + sum(DEV.NPML);
Sy      = Ny * dy;

% COMPUTE 2X GRID
Nx2     = 2 * Nx;
dx2     = dx/2;
Ny2     = 2 * Ny;
dy2     = dy/2;
DEV.RES = [dx2 dy2];

% COMPUTE AXES
xa      = [0:Nx-1]*dx;
xa      = xa - mean(xa);
ya      = [0:Ny-1]*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON 2X GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATERIAL ARRAYS
DEV.ER2     = er1 * ones(Nx2,Ny2);
DEV.UR2     = ur1 * ones(Nx2,Ny2);

% COMPUTE START AND STOP INDICES
% X Indices
nx  = 1;
nx1 = ceil(Nx2/2) - 1;
nx2 = nx1  - 1;
nx3 = Nx2;

% Y INDICES
ny  = ceil(Ny2/2);
ny1 = ny + ceil(T/dx2) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPLEMENT FDFD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE AUXILIAY ARRAYS
REF = zeros(1,NLAM);
TRN = zeros(1,NLAM);
CON = zeros(1,NLAM);

display('============================');
display('Now Entering FDFD Simulation');
display('============================');

% PARAMETER SWEEP LOOP
for i = 1 : NLAM
    
    tic
    % CHANGE VALUE OF FREQUENCY INSIDE SOURCE STRUCTURE
    SRC.lam0 = lam0(i);
    
    % PERFORM FDFD ANALYSIS
    DAT = fdfd2d(DEV,SRC);
    clc;
    display('============================================');
    display(['Finished ' num2str(i) 'th Backwards Divide']);
    display('============================================');
    REF(i) = DAT.REF;
    TRN(i) = DAT.TRN;
    CON(i) = DAT.CON;
    time = toc;
    disp(['Estimated Time: ' num2str(time*(NLAM-i)/60) ' minutes']);
    disp([num2str(i) ' out of ' num2str(NLAM) ' Iterations']);
    
end

% DISPLAY FINAL ANSWER
a = plot(lam0,REF,'LineWidth',1.5,'Color','r');
hold on;
a2 = plot(lam0,TRN,'LineWidth',1.5,'Color','b');
a3 = plot(lam0,CON,'LineWidth',1.5,'Color','k');
hold off;
a4 = get(a(1),'Parent');
set(a4,'FontSize',11);
title([ MODE ' Mode Simulation at ' num2str(theta) '$^{\circ}$'],...
    'Interpreter','LaTex','FontSize',16);
xlabel('$\textrm{Frequency (GHz)}$','Interpreter','LaTex','FontSize',12);
ylabel('$\textrm{Amplitude}$','Rotation',90,'Interpreter',...
    'LaTex','FontSize',12);
ylim([0 105]);
a = legend('Reflectance','Transmittance','Conservation');
set(a,'Box','off','FontSize',10.5,'Location','West');
