% Optimization to find effective permittivity for grating layer
clc;
clear all;
close all;

% Units
meters = 1;
micrometers = 1e-6 * meters;
seconds = 1;

c0 = 299792458 * meters/seconds;       % Speed of light
u0 = 4*pi * 1e-7 * 1/meters;           % Free space permeability
e0 = 1/(u0*c0^2) * 1/meters;       % Free space permittivity
N0 = sqrt(u0/e0);                      % Free space impedance

% Equations
ref  =@(N2,N1) (N2-N1)./(N2+N1);
cref =@(r12,r23,psi) (r12+r23.*exp(-2i*psi))./(1+r12.*r23.*exp(-2i.*psi));
N    =@(ur,er) sqrt(ur*u0/er/e0);

% Dashboard for TMM
lam0 = 1.55 * micrometers;
k0   = 2*pi/lam0;

% Device parameters
DEV.er1 = 1.0;
DEV.ur1 = 1.0;
DEV.er2 = 1.0;
DEV.ur2 = 1.0;

DEV.ER = [1.5 2.0];
DEV.UR = [1.0 1.0];
DEV.L  = [0.001 1.0]*lam0;
DEV.NP = 2;

% Source parameters
SRC.lam0 = lam0;
SRC.theta = 0;
SRC.phi = 0;
SRC.ate = 1;
SRC.atm = 0;

% Optimize
tol = 1e-6;
err = inf;
limit = 1000; lim = 0;
gdstep = 0.01 * lam0;
gestep = 0.01 * DEV.ER(1);
alpha = 0.01;

while (err > tol)
    
    % Calculate reflectance
    DAT = tmm1d(DEV,SRC);
    disp(['Reflectance: ' num2str(DAT.REF*100) ' %']);
    
    % Calulate gradient around value
    grd = [DEV.L(1)+gdstep  , DEV.L(1)+gdstep , DEV.L(1)+gdstep;
           DEV.L(1)         , DEV.L(1)        , DEV.L(1)       ;
           DEV.L(1)-gdstep  , DEV.L(1)-gdstep , DEV.L(1)-gdstep];
       
    gre = [DEV.ER(1)-gestep , DEV.ER(1)       , DEV.ER(1)+gestep;
           DEV.ER(1)-gestep , DEV.ER(1)       , DEV.ER(1)+gestep;
           DEV.ER(1)-gestep , DEV.ER(1)       , DEV.ER(1)+gestep];
    
    DEVT = DEV;
    f = ones(size(grd));
    
    for n = 1:length(grd(:))
        DEVT.L(1)  = grd(n);
        DEVT.ER(1) = gre(n);
        DATT = tmm1d(DEVT,SRC);
        f(n) = DATT.REF;
    end
    
    grad = f - DAT.REF;
    
    % Find maximum value index
    ind = find(~(grad-max(max(grad))));

    % Move in the direction of the gradient
    DEV.L(1)  = DEV.L(1)  + alpha*grd(ind);
    DEV.ER(1) = DEV.ER(1) + alpha*gre(ind);
    
    % Calculate error
    err = abs(DAT.REF - DATT.REF)
    lim = lim + 1;
    
    if lim == limit
        break;
    end
end

disp(['Grating depth = ' num2str(DEV.L(1)/lam0) ' wavelengths']);
disp(['Grating permittivity = ' num2str(DEV.ER(1))]);

















