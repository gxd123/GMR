function [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc)
% YEEDER Construct Yee Grid Derivative Operators on a 2D Grid
%
% [DEX,DEY,DHX,DHY] = yeeder(NGRID,RES,BC,kinc);
%
% Note for normalized grid, use this function as follows:
%
% [DEX,DEY,DHX,DHY] = yeeder(NGRID,k0*RES,BC,kinc/k0);
%
% Input Arguments
% =================
% NGRID [Nx Ny] grid size
% RES [dx dy] grid resolution of the 1X grid
% BC [xbc ybc] boundary conditions
% -2: periodic (requires kinc)
% 0: Dirichlet
% kinc [kx ky] incident wave vector
% This argument is only needed for periodic boundaries.

A = NGRID(1)*NGRID(2);  % Size of Derivative Matrices

% CHECK IF BOUNDARY CONDITIONS AND KINC ARE DEFINED
if(~exist('BC','var'))
    BC = [0 0];
end
if(~exist('kinc','var'))
    kinc = [0 0];
end

% COMPUTE BASED ON GRID SIZES
if NGRID(1) == 1
    if(~BC(1))
        DEX = sparse(A,A);
    elseif(BC(1) == -2)
        DEX = 1i*kinc(1)*speye(A);
    end
    DHX = DEX;
else
    % Do middle diagonal first
    DEX = -speye(A);
    
    % Now the second diagonal based on boundary conditions
    d = ones(A,1);
    for n = 1:NGRID(1):A
        d(n) = 0;
    end
    DEX = spdiags(d,1,DEX);
    if(BC(1) == -2)
        d = zeros(A,1);
        for n = 1:NGRID(1):A
            d(n) = exp(1i*kinc(1)*(NGRID(1)*RES(1)));
        end
        DEX = spdiags(d,-(NGRID(1)-1),DEX);
    end
    DEX = DEX/RES(1);
    DHX = -DEX';
end

if NGRID(2) == 1
    if(~BC(2))
        DEY = sparse(A,A);
    elseif(BC(2) == -2)
        DEY = 1i*kinc(2)*speye(A);
    end
    DHY = DEY;
else
    % Do middle diagonal first
    DEY = -speye(A);
    
    % Now the second diagonal based on boundary conditions
    d = ones(A,1);
    DEY = spdiags(d,NGRID(1),DEY);
    if(BC(2) == -2)
        d = exp(1i*kinc(2)*(NGRID(2)*RES(2)))*ones(A,1);
        DEY = spdiags(d,-(A-NGRID(1)),DEY);
    end
    DEY = DEY/RES(2);
    DHY = -DEY';
end