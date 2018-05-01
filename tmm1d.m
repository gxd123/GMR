function DAT = tmm1d(DEV,SRC)
% TMM1D One-Dimensional Transfer Matrix Method
%
% DAT = tmm1d(DEV,SRC);
%
% Homework #5, Problem #3
% EE 5337 - COMPUTATIONADEV.L EDEV.LECTROMAGNETICS
%
% IDEV.NPUT ARGUMENTS
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

% UNITS
j = 1i;
I = [1 0;0 1];
O = [0 0;0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SODEV.URCE PARAMETDEV.ERS
SRC.theta = SRC.theta;          % Elevation angle in degrees
SRC.phi   = SRC.phi;            % Azimuthal angle in degrees

% DEFINE LAYERS
% DEV.L  = DEV.L*SRC.lam0;                  % Array of the thickness of each layer
                                          % relative to wavelength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPLEMENT TRANSFER MATRIX METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute wave vector components
k0 = 2*pi/SRC.lam0; 
kx = sqrt(DEV.ur1*DEV.er1)*sin(SRC.theta)*cos(SRC.phi);         % Normalized x component
ky = sqrt(DEV.ur1*DEV.er1)*sin(SRC.theta)*sin(SRC.phi);         % Normalized y component

% Compute gap medium parameters 
Qh = [kx*ky 1+ky^2 ; -(1+kx^2) -kx*ky];
Vh = -j*Qh;

% Initialize global scattering matrix
SG = struct('S11',O,'S12',I,'S21',I,'S22',O);
Si = struct('S11',[],'S12',[],'S21',[],'S22',[]);

% Main loop for each layer of the device
for n = 1:length(DEV.UR)
    % Calculate Parameters for Layer i 
    Q = [kx*ky DEV.UR(n)*DEV.ER(n)-kx^2;ky^2-DEV.UR(n)*DEV.ER(n) -kx*ky]/DEV.UR(n);
    kz = sqrt(DEV.UR(n).*DEV.ER(n) - kx^2 - ky^2);       % Normalized kz vector
    OMEGA = j*kz*I;
    V = Q/OMEGA;
    
    % Calculate Scattering Matrix for Layer i
    X = diag(exp(diag(OMEGA)*k0*DEV.L(n)));
    A = I+V\Vh;
    B = I-V\Vh;
    D = A-X*B*(A\X)*B;
    Si.S11 = D\(X*B*(A\X)*A-B);
    Si.S12 = D\X*(A-B*(A\B));
    Si.S21 = Si.S12;
    Si.S22 = Si.S11;
    
    % Update Global Scattering Matrix
    SG = star(SG,Si);
end

% Cascade for repeated cells
SG = cascn(SG,DEV.NP);

% Connect to external regions
% First reflection region (don't use star because of order)
Qref = [kx*ky DEV.ur1*DEV.er1-kx^2;ky^2-DEV.ur1*DEV.er1 -kx*ky]/DEV.ur1;
kzref = sqrt(DEV.ur1.*DEV.er1 - kx^2 - ky^2);       % Normalized kz reflection region vector
OMEGA = j*kzref*I;
Vref = Qref/OMEGA;
Aref = I+Vh\Vref;
Bref = I-Vh\Vref;

SR11 = -Aref\Bref;
SR12 = 2*inv(Aref);
SR21 = 0.5*(Aref-Bref/Aref*Bref);
SR22 = Bref/Aref;

% Now transmission region
Qtrn = [kx*ky DEV.ur2*DEV.er2-kx^2;ky^2-DEV.ur2*DEV.er2 -kx*ky]/DEV.ur2;
kztrn = sqrt(DEV.ur2.*DEV.er2 - kx^2 - ky^2);       % Normalized kz transmission region vector
OMEGA = j*kztrn*I;
Vtrn = Qtrn/OMEGA;
Atrn = I+Vh\Vtrn;
Btrn = I-Vh\Vtrn;

ST11 = Btrn/Atrn;
ST12 = 0.5*(Atrn-Btrn/Atrn*Btrn);
ST21 = 2*inv(Atrn);
ST22 = -Atrn\Btrn;

% Update global scattering matrix
% First reflection region
D = SR12/(I-SG.S11*SR22);
F = SG.S21/(I-SR22*SG.S11);
SG.S22 = SG.S22 + F*SR22*SG.S12;
SG.S21 = F*SR21;
SG.S12 = D*SG.S12;
SG.S11 = SR11 + D*SG.S11*SR21;

% Now transmission region
D = SG.S12/(I-ST11*SG.S22);
F = ST21/(I-SG.S22*ST11);
SG.S11 = SG.S11 + D*ST11*SG.S21;
SG.S12 = D*ST12;
SG.S21 = F*SG.S21;
SG.S22 = ST22 + F*SG.S22*ST12;

% Calculate the source
kinc = k0*[kx;ky;sqrt(DEV.ur1*DEV.er1)*cos(SRC.theta)];
n = [0;0;1];
if(SRC.theta ~= 0)
    aTE = cross(n,kinc)/norm(cross(n,kinc));
else
    aTE = [0;1;0];
end
aTM = cross(aTE,kinc)/norm(cross(aTE,kinc));
P = SRC.ate*aTE + SRC.atm*aTM;
Esrc = [P(1);P(2)]./norm(P);

% Calculate transmitted and refelcted fields
Eref = SG.S11*Esrc;
Etrn = SG.S21*Esrc;

% Calculate longitudinal field components
Ezref = -(kx*Eref(1)+ky*Eref(2))/kzref;
Eztrn = -(kx*Etrn(1)+ky*Etrn(2))/kztrn;

% Calculate Tranmittance and Reflectance
Eref = [Eref;Ezref];
Etrn = [Etrn;Eztrn];
DAT.REF = norm(Eref)^2;
DAT.TRN = norm(Etrn)^2*(real(k0*kztrn/DEV.ur2))/(real(kinc(3)/DEV.ur1));
DAT.CON = DAT.REF + DAT.TRN;

end