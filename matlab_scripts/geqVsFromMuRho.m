function Vs = geqVsFromMuRho(Mu, Rho)
%% function Vs = geqVsFromMuRho(Mu, Rho)
% calculate S-wave velocity: Vs = (Mu./Rho).^0.5.
% 
%   Vs              - S-wave velocity [km/s]
%   G            	- shear modulus [GPa]
%   rho             - density  [g/ccm]
% 
%_______________________________________________________________________
% Examples
% Vs    = geqVsFromMuRho(30, 2.65)
% Vs    = geqVsFromMuRho([30 35], 2.65)
% Vs    = geqVsFromMuRho([30 35], [2.65 2.67])
%_______________________________________________________________________
% Erling Hugo Jensen, 08/09/09
%
% See also geqMuFromRhoVpVs, geqVpFromKMuRho

%% Check input parameters
narginchk(2,2);

%% Calculate 
Vs = (Mu./Rho).^0.5;
