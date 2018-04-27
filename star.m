function S = star(SA,SB)
% STAR Redheffer Star Product
%
% S = star(SA,SB)
%
% INPUT ARGUMENTS
% ================
% SA        First Scattering Matrix
%   .S11    S11 scattering parameter
%   .S12    S12 scattering parameter
%   .S21    S21 scattering parameter
%   .S22    S22 scattering parameter
%
% SB        Second Scattering Matrix
%   .S11    S11 scattering parameter
%   .S12    S12 scattering parameter
%   .S21    S21 scattering parameter
%   .S22    S22 scattering parameter
%
% OUTPUT ARGUMENTS
% ================
% S         Combined Scattering Matrix
%   .S11    S11 scattering parameter
%   .S12    S12 scattering parameter
%   .S21    S21 scattering parameter
%   .S22    S22 scattering parameter

S  = struct('S11',[],'S12',[],'S21',[],'S22',[]);
% SA = struct('S11',SA(1:length(SA)/2,1:length(SA)/2),...
%             'S12',SA(1:length(SA)/2,1+length(SA)/2:end),...
%             'S21',SA(1+length(SA)/2:end,1:length(SA)/2),...
%             'S22',SA(1+length(SA)/2:end,1+length(SA)/2:end));
% SB = struct('S11',SB(1:length(SB)/2,1:length(SB)/2),...
%             'S12',SB(1:length(SB)/2,1+length(SB)/2:end),...
%             'S21',SB(1+length(SB)/2:end,1:length(SB)/2),...
%             'S22',SB(1+length(SB)/2:end,1+length(SA)/2:end));
I = eye(length(SA.S11));
D = SA.S12/(I-SB.S11*SA.S22);
F = SB.S21/(I-SA.S22*SB.S11);
S.S11 = SA.S11 + D*SB.S11*SA.S21;
S.S12 = D*SB.S12;
S.S21 = F*SA.S21;
S.S22 = SB.S22 + F*SA.S22*SB.S12;
end