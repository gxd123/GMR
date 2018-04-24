function DAT = fdfd2d(DEV,SRC)
% FDFD2D Two-Dimensional Finite-Difference Frequency-Domain
%
% DAT = fdfd2d(DEV,SRC)
%
% INPUT ARGUMENTS
% =================
% DEV Device Parameters
% .UR2 Relative permeability on 2X grid
% .ER2 Relative permittivity on 2X grid
% .NPML Size of PML on 1X grid [xlo xhi ylo yhi]
% .RES [dx2 dy2] grid resolution of 2X grid
%
% SRC
% .lam0 free space wavelength
% .theta Angle of incidence
% .MODE Mode: 'E' or 'H'
%
% OUTPUT ARGUMENTS
% ================= 
% DAT Output Data
% .RDE Array of diffraction efficiencies of reflected harmonics
% .REF Overall Reflectance
% .TDE Array of diffraction efficiencies of transmitted harmonics
% .TRN Overall Transmittance
% .CON Conservation of Energy
% .F Field
%
% Homework #9, Problem 1
% EE 5320 - COMPUTATIONAL ELECTROMAGNETICS

% Get grid sizes and resolutions
if(size(DEV.UR2) ~= size(DEV.ER2))
    error('UR2 and ER2 are not the same size');
end
[Nx2,Ny2] = size(DEV.UR2);
Nx  = Nx2/2; Ny = Ny2/2;
dx  = DEV.RES(1)*2;
dy  = DEV.RES(2)*2;
ur1 = DEV.UR2(1,1); ur2 = DEV.UR2(end,end);
er1 = DEV.ER2(1,1); er2 = DEV.ER2(end,end);

% Get PML and incorporate it
[sx,sy] = calcpml2d([Nx2 Ny2],2*DEV.NPML);
URxx = DEV.UR2./sx.*sy;
URyy = DEV.UR2.*sx./sy;
URzz = DEV.UR2.*sx.*sy;
ERxx = DEV.ER2./sx.*sy;
ERyy = DEV.ER2.*sx./sy;
ERzz = DEV.ER2.*sx.*sy;

% Overlay materials onto 1x grid
URxx = URxx(1:2:Nx2,2:2:Ny2);
URyy = URyy(2:2:Nx2,1:2:Ny2);
URzz = URzz(2:2:Nx2,2:2:Ny2);
ERxx = ERxx(2:2:Nx2,1:2:Ny2);
ERyy = ERyy(1:2:Nx2,2:2:Ny2);
ERzz = ERzz(1:2:Nx2,1:2:Ny2);

% Compute the wave vector terms
nref = sqrt(ur1*er1);
ntrn = sqrt(ur2*er2);
m = [-floor(Nx/2):floor(Nx/2)]';
k0 = 2*pi/SRC.lam0;
kinc = k0*nref*[sin(SRC.theta) cos(SRC.theta)];
kx = kinc(1) - m*2*pi/(dx*Nx);
kyref = sqrt((k0*nref).^2 - kx.^2);
kytrn = sqrt((k0*ntrn).^2 - kx.^2);

% Construct diagonal materials matrices
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));

% Construct derivative matrices
[DEX,DEY,DHX,DHY] = yeeder([Nx Ny],k0*[dx dy],[-2 0],kinc/k0);

% Compute the wave matrix A
if(SRC.MODE == 'E')
    A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
elseif(SRC.MODE == 'H');
    A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
else
    error('Invalid Mode');
end

% Compute the source field
xa = dx.*m;
ya = dy.*[0:Ny-1];
[Y,X] = meshgrid(ya,xa);
fsrc = exp(1i*(kinc(1)*X + kinc(2)*Y));
fsrc = fsrc(:);

% Compute scattered-field masking matrix
Q = zeros(Nx,Ny);
Q(:,1:(DEV.NPML(3)+100)) = 1;
Q = diag(sparse(Q(:)));

% Compute the source vector
b = (Q*A - A*Q)*fsrc;

% Compute the field
f = A\b;
f = reshape(f,Nx,Ny);
DAT.F = f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Transmitted and Reflected Fields
Eref = f(:,DEV.NPML(3)+1);
Etrn = f(:,end-DEV.NPML(4));

% Remove phase tilt
Aref = Eref.*exp(-1i*kinc(1).*xa);
Atrn = Etrn.*exp(-1i*kinc(1).*xa);

% Calculate complex amplitudes of the spatial harmonics
Sref = flipud(fftshift(fft(Aref)))/Nx;
Strn = flipud(fftshift(fft(Atrn)))/Nx;

% Calculate diffraction efficiencies (for source wave amplitude of 1)
if(SRC.MODE == 'E')
    DAT.RDE = abs(Sref).^2 .* real(kyref./ur1)./real(kinc(2)./ur1);
    DAT.TDE = abs(Strn).^2 .* real(kytrn./ur2)./real(kinc(2)./ur1);
elseif(SRC.MODE == 'H')
    DAT.RDE = abs(Sref).^2 .* real(kyref./er1)./real(kinc(2)./er1);
    DAT.TDE = abs(Strn).^2 .* real(kytrn./er2)./real(kinc(2)./er1);
end

% Calculate Reflectance, Transmittance, and Conservation of Power
DAT.REF = sum(DAT.RDE);
DAT.TRN = sum(DAT.TDE);
DAT.CON = DAT.REF + DAT.TRN;


