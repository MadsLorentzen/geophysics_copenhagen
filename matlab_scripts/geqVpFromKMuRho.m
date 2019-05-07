function Vp = geqVpFromKMuRho(K, Mu, Rho)
%% function Vp = geqVpFromKMuRho(K, Mu, Rho)
% calculate P-wave velocity: Vp = ((K + 4/3*.Mu)./Rho).^(0.5).
% 
%   Vp              - P-wave velocity [km/s]
%   K            	- bulk modulus [GPa]
%   G               - shear modulus [GPa]
%   rho             - density [g/ccm]
% 
%_______________________________________________________________________
% Examples
% Vp    = geqVpFromKMuRho(36, 30, 2.65)
% Vp    = geqVpFromKMuRho([36 38], 30, 2.65)
% Vp    = geqVpFromKMuRho([36 38], [30 33], 2.65)
% Vp    = geqVpFromKMuRho([36 38], [30 33], [2.65 2.67])
%_______________________________________________________________________
% Erling Hugo Jensen, 08/09/09
%
% See also geqKFromRhoVpVs, geqMuFromRhoVpVs, geqVsFromKMuRho.

%% Check input parameters
narginchk(3,3);

%% Calculate 
Vp = ((K + 4/3*Mu)./Rho).^0.5;
