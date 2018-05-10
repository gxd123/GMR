function SC = cascn(S,N)
% CASC Cascading and Doubling Algorithm
%
% SC = cascn(S,N)
%
% This MATLAB function uses an efficient doubling algorithm
% to cascade N periods.
%
% INPUT ARGUMENTS
% ================
% S         Scattering Matrix for One Period
%   .S11    S11 scattering parameter
%   .S12    S12 scattering parameter
%   .S21    S21 scattering parameter
%   .S22    S22 scattering parameter
%
% N         Number of scattering matrices to cascade
%
% OUTPUT ARGUMENTS
% ================
% SC        Overall Scattering Matrix for Cascade
%   .S11    S11 scattering parameter
%   .S12    S12 scattering parameter
%   .S21    S21 scattering parameter
%   .S22    S22 scattering parameter

% Create common I and O matrices
I = eye(length(S.S11));
O = zeros(length(S.S11),length(S.S11));

% Initialize matrices
SC  = struct('S11',O,'S12',I,'S21',I,'S22',O);
Sb = S;

% Loop through cascading
n = de2bi(N);   % Convert to binary array
for k = 1:length(n)
    if(n(k)==1)
        SC = star(SC,Sb);   % If 1, update output
    end
    if(k ~= length(n))
        Sb = star(Sb,Sb);       % Double Sbin        
    end
end

