%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMR_Resonance.m
% This script simulates a 1D GMR Filter to find its resonant wavelength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','docked');

% RESTORE STATE
clear all;
close all;
clc;

% OPEN FIGURE WINDOW
fig  = figure('Color','w');

% UNITS
micrometers = 1;
millimeters = 1e3 * micrometers;
centimeters = 1e4 * micrometers;
meters      = 1e3 * centimeters;
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
NFREQ       = 501;
NLAM        = NFREQ;
lam1        = 1.4 * micrometers;
lam2        = 1.7 * micrometers;
lam0        = linspace(lam1,lam2,NLAM);
SRC.theta   = 0 * degrees;
SRC.MODE    = 'E';

% GRATING PARAMETERS
lamd    = 1.55 * micrometers;
fd      = c0/lamd;
er      = 2.0;
ur      = 1.0;
nr      = sqrt(ur*er);
t       = lamd/(2*nr);              % Substrate thickness 
L       = 0.856092903225807 * lamd; % Grating period
% L       = 0.6 * lamd;
d       = 0.152547096774194 * lamd; % Grating depth
% d       = 0.15 * lamd;
f       = 0.454096;                     % Duty cycle
% f       = 0.8;

% EXTERNAL MATERIALS
ur1 = 1.0; %permeability in the reflection region
er1 = 1.0; %permittivity in the reflection region
ur2 = 1.0; %permeability in the transmission region
er2 = 2.25; %permittivity in the transmission region

% GRID PARAMETERS
NRES = 100;
BUFZ = 2*max(lam0) * [1 1];
DEV.NPML = [0 0 20 20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE REFRACTIVE INDICES
nr1 = sqrt(ur1 * er1);
nr2 = sqrt(ur2 * er2);
nr  = sqrt(ur * er);
nmax    = max([nr nr1 nr2]);

% RESOLVE FOR MINIMUM WAVELENGTH
lam_min = min(lam0);

% COMPUTE INITIAL GRID RESOLUTION
dx      = lam_min/nmax/NRES;
dy      = lam_min/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSION
Nx      = 2*ceil(L/dx/2) + 1; % Ensure grid is always odd
dx      = L/Nx;
Ny      = ceil((d+t)/dy); 
dy      = (d+t)/Ny; % Recalculate resolution

% COMPUTE GRID SIZE
Sx      = L;
Sy      = d + t + sum(BUFZ);
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
wx  = round(f*Nx2);
nx1 = round((Nx2-wx)/2);
nx2 = nx1 + wx - 1;

% Y Indices
ny1 = 2*DEV.NPML(3) + round(BUFZ(1)/dy2) + 1;
ny2 = ny1 + round(d/dy2) - 1;
ny3 = ny2 + round(t/dy2) - 1;
ny4 = ny2 + round(d/dy2) - 1;

% INCORPORATE GRATING
DEV.ER2(:,ny2:ny3)          = er;
DEV.ER2(nx1:nx2,ny2:ny4)    = er1;
DEV.ER2(:,ny3:Ny2)          = er2;

imagesc(DEV.ER2');
axis equal tight;
colormap('Gray');
colorbar;

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
    
    subplot(121);
    plot(lam0,REF);
    
    subplot(122)
    imagesc(xa,ya,DEV.ER2');
    hold on
    imagesc(xa,ya,abs(DAT.F'));
    hold off;
    axis equal tight
    colormap('Jet');
    drawnow
    
end

% DISPLAY FINAL ANSWER
figure('Color','w');
a = plot(lam0,REF,'LineWidth',1.5,'Color','r');
hold on;
a2 = plot(lam0,TRN,'LineWidth',1.5,'Color','b');
a3 = plot(lam0,CON,'LineWidth',1.5,'Color','k');
hold off;
a4 = get(a(1),'Parent');
set(a4,'FontSize',11);
title([ SRC.MODE ' Mode Simulation at ' num2str(SRC.theta) '$^{\circ}$'],...
    'Interpreter','LaTex','FontSize',16);
xlabel('$\textrm{Wavelength $\lambda_0$}$','Interpreter','LaTex','FontSize',12);
ylabel('$\textrm{Amplitude}$','Rotation',90,'Interpreter',...
    'LaTex','FontSize',12);
ylim([0 1.05]);
xlim([lam1 lam2]);
a = legend('Reflectance','Transmittance','Conservation');
set(a,'Box','off','FontSize',10.5,'Location','West');

% COMPUTE FWHM
MAX_REF     = max(REF);
MAX_IND     = find(REF == MAX_REF);
H_MAX       = 0.5 * MAX_REF;
HM_INDEX    = find(REF > H_MAX);
FWHM        = lam0(max(HM_INDEX)) - lam0(min(HM_INDEX));