%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMRFFSweep.m
% 
% This MATLAB script does an extensive search to maximize reflectance of a
% GMR filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','docked');

% RESTORE STATE
clc;
clear all;
close all;

% UNITS
micrometers = 1;
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
NFREQ       = 501;
SRC.lam0    = 1.55 * micrometers;    % Free space wavelength
SRC.theta   = 0 * degrees;
SRC.MODE    = 'E';                   % EM mode
f0          = c0/SRC.lam0;
lam1       = 1.4 * micrometers;
lam2       = 1.7 * micrometers;
lam0        = linspace(lam1,lam2,NFREQ);

% GRATING PARAMETERS
lamd = 1.55 * micrometers;       % Design wavelength
fd   = c0/lamd;                  % Design frequency
ur   = 1.0;                      % Grating permeability
er   = 2.0;                      % Grating permittivity
nr   = sqrt(ur*er);              % Substrate refractive index 
L    = 0.856092903225807 * lamd; % Grating period
d    = 0.152547096774194 * lamd; % Grating depth
t    = lamd/(2*nr);              % Substrate thickness

% EXTERNAL MATERIALS
ur1 = 1.0;                    % Reflection region permeability
er1 = 1.0;                    % Reflection region permittivity
ur2 = 1.0;                    % Transmission region permeability
er2 = 1.0;                    % Transmission region permittivity

% GRID PARAMETERS
NRES = 100;                    % Grid resolution
BUFZ = 2*lam2 * [1 1];       % Spacer region above and below grating
DEV.NPML = [ 0 0 20 20];           % Size of PML at top and bottom of grid

% INITIALIZE ARRAYS
REF = zeros(1,NFREQ);
TRN = REF;
CON = REF;

% PERIOD SWEEP
fr = linspace(0.001,0.9,NFREQ);

VIS = 1;

for n = 1:NFREQ
    f    = fr(n);             % Grating period
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CALCULATE OPTIMIZED GRID
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BUILD DEVICE ON 2x GRID
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DEV.UR2 = ur1*ones(Nx2,Ny2);
    DEV.ER2 = er1*ones(Nx2,Ny2);
    
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
    DEV.ER2(nx1:nx2,ny2:ny4)    = 1.0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% IMPLEMENT FDFD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    SRC.lam0 = lamd;         %angle of incidence
    DAT = fdfd2d(DEV,SRC);
    REF(n) = DAT.REF;
    TRN(n) = DAT.TRN;
    CON(n) = DAT.CON;
    clc;
    time = toc;
    mint = floor(time*(NFREQ-n)/60);
    sec = round((time*(NFREQ-n)/60 - mint)*60);
    disp(['Estimated Time: ' num2str(mint) ' minutes '...
        num2str(sec) ' seconds']);
    disp([num2str(n) ' out of ' num2str(NFREQ) ' Iterations']);
    
    if VIS
        plot(fr,REF);
        title([SRC.MODE ' Mode @ ' num2str(SRC.lam0) ' \mum']);
        shading interp;
        xlabel('Fill Fraction'); ylabel('REF');
        drawnow;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(fr,REF,'LineWidth',1.5);
title('Period Sweep','FontSize',14);
ylabel('Reflectance','FontSize',12);
xlabel('$\frac{L}{\lambda_0}$','FontSize',12,'Interpreter','LaTex');

% OBTAIN BEST FIT FOR PERIOD
A = (max(unique(REF)));
A = find(REF == A);
f = fr(A);
