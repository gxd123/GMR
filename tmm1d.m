function DAT = tmm1d(DEV,SRC)
% TMM1D One-Dimensional Transfer Matrix Method
%
% DAT = tmm1d(DEV,SRC);
%
% Homework #5, Problem #3
% EE 5337 - COMPUTATIONAL ELECTROMAGNETICS
%
% INPUT ARGUMENTS
% ================
% DEV Device Parameters
%
% .er1 relative permittivity in reflection region
% .ur1 relative permeability in reflection region
% .er2 relative permittivity in transmission region
% .ur2 relative permeability in transmission region
%
% .ER array containing permittivity of each layer in unit cell
% .UR array containing permeability of each layer in unit cell
% .L array containing thickness of each layer in unit cell
% .NP number of unit cells to cascade
%
% SRC Source Parameters
%
% .lam0 free space wavelength
%
% .theta elevation angle of incidence (radians)
% .phi azimuthal angle of incidence (radians)
%
% .ate amplitude of TE polarization
% .atm amplitude of TM polarization
%
% OUTPUT ARGUMENTS
% ================
% DAT Output Data
%
% .REF Reflectance
% .TRN Transmittance

% MISC HOUSEKEEPING
I       = eye(2,2);
Z       = zeros(2,2);

% DEFINE SIMULATION PATAMETERS
theta   = SRC.theta;
phi     = SRC.phi;
pte     = SRC.pte;
ptm     = SRC.ptm;
ur1     = DEV.ur1;
er1     = DEV.er1;
ur2     = DEV.ur2;
er2     = DEV.er2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPLEMENT TRANSFER MATRIX METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Wave Vector Components
k0      = 2*pi/SRC.lam0;
ninc    = sqrt(ur1*er1);
kx      = ninc*sin(theta)*cos(phi);
ky      = ninc*sin(theta)*sin(phi);

% Compute medium gap parameters
Qh = [ kx*ky 1+kx^2 ; -(1+kx^2) -kx*ky ];
Vh = -1i*Qh;

% Initialize global s matrix
SG = struct('S11',Z,'S12',I,'S21',I,'S22',Z);
Si = struct('S11',[],'S12',[],'S21',[],'S22',[]);

%
% Main Loop
%
NLAY = length(DEV.L);
for nlay = 1 : NLAY
    
        % CALCULATE PARAMETERS FOR ith LAYER
        ur  = DEV.UR(nlay);
        er  = DEV.ER(nlay);
        l   = DEV.L(nlay);
        
        % SORT EIGEN MODES FOR nlay^th layer
        kz  = sqrt(ur*er - kx^2 -ky^2);
        Q   = (1/ur)*[ kx*ky    ur*er-kx^2 ; ky^2-ur*er   -kx*ky];
        OMEGA = 1i*kz*I;
        V   = Q/OMEGA;
        
        % CALCULATE SCATERING MATRIX FOR LAYER
        A = I + V\Vh;
        B = I - V\Vh;
        X = diag(exp(diag(OMEGA)*k0*l));
        D = A - (X*B/A)*X*B;
        Si.S11  = D\(((X*B/A)*X*A)-B);
        Si.S12  = D\X*(A-(B/A)*B);
        Si.S21  = Si.S12;
        Si.S22  = Si.S11;
        
        % Update global scattering matrix
        SG = star(SG,Si);
        
end

% CASCADE FOR REPEATED CELLS
SG = cascn(SG,DEV.NP);

% COMPUTE EIGEN MODES IN REFLECTION REGION
Q = (1/ur1)*[kx*ky ur1*er1-kx^2; ky^2-ur1*er1 -kx*ky];
kzref = sqrt(ur1*er1 - kx^2 - ky^2);
OMEGA = 1i*kzref*I;
Vref = Q/OMEGA;

% CALCULATE SCATTERING MATRIX TO CONNECT TO REFLECTION REGION
A = I + Vh\Vref;
B = I - Vh\Vref;
SR11 = -A\B;
SR12 = 2*I/A;
SR21 = 0.5*(A - B/A*B);
SR22 = B/A;

% COMPUTE EIGEN MODES IN TRANSMISSION REGION
Q = (1/ur2)*[kx*ky ur2*er2-kx^2; ky^2-ur2*er2 -kx*ky];
kztrn = sqrt(ur2*er2 - kx^2 - ky^2);
OMEGA = 1i*kztrn*I;
Vtrn = Q/OMEGA;

% CALCULATE TRANSMISSION SIDE SCATTERING MATRIX
A = I + Vh\Vtrn;
B = I - Vh\Vtrn;
ST11 = B/A;
ST12 = 0.5*(A - B/A*B);
ST21 = 2*I/A;
ST22 = -A\B;

% UPDATE GLOBAL SCATTERING MATRIX
D = SR12/(I-SG.S11*SR22);
F = SG.S21/(I-SR22*SG.S11);
SG.S22 = SG.S22 + F*SR22*SG.S12;
SG.S21 = F*SR21;
SG.S12 = D*SG.S12;
SG.S11 = SR11 + D*SG.S11*SR21;

D = SG.S12/(I-ST11*SG.S22);
F = ST21/(I-SG.S22*ST11);
SG.S11 = SG.S11 + D*ST11*SG.S21;
SG.S12 = D*ST12;
SG.S21 = F*SG.S21;
SG.S22 = ST22 + F*SG.S22*ST12;

% CALCULATE THE SOURCE
kinc = k0*ninc*[sin(theta)*cos(phi) ; sin(theta)*sin(phi) ; cos(theta)];
% Surface Normal
n = [0;0;1];
if theta == 0
    ate = [0;1;0];
else
    ate = cross(n,kinc);
    ate = ate/norm(ate);
end
atm = cross(ate,kinc);
atm = atm/norm(atm);
P = pte*ate + ptm*atm;
P = P/norm(P);

% CALCULATE COEFFICIENTS FOR INCIDENT FIELD
Esrc = [ P(1) ; P(2) ];
csrc = I\Esrc;

% COMPUTE REFLECTED AND TRANSMITTED FIELDS
Eref = I*SG.S11*csrc;
Etrn = I*SG.S21*csrc;

% COMPUTE LONGITUDINAL FIELD COMPONENTS
Exref = Eref(1);
Eyref = Eref(2);
Ezref = - (kx*Exref + ky*Eyref)/kzref;
Eref = [ Exref ; Eyref ; Ezref ];

Extrn = Etrn(1);
Eytrn = Etrn(2);
Eztrn = - (kx*Extrn + ky*Eytrn)/kztrn;
Etrn = [ Extrn ; Eytrn ; Eztrn ];

% CALCULATE REFLECTANCE AND TRANSMITTANCE
DAT.REF = (norm(Eref)^2);
DAT.TRN = norm(Etrn)^2 * real(ur1/ur2*kztrn/kzref);