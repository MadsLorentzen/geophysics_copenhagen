function nu = geqPoissonFromKMu(K, Mu)
%% function poissonsRatio = geqPoissonFromKMu(K, Mu)
% calculate poisson's ratio: nu = {3*K - 2*Mu}./{2*(3*K + Mu)}.
% 
%   nu              - Poisson's ratio
%   K               - Bulk modulus  
%   Mu              - Shear modulus 
% 
%_______________________________________________________________________
% Examples
% nu    = geqPoissonFromKMu(15, 21)
% nu    = geqPoissonFromKMu([15 19], 21)
% nu    = geqPoissonFromKMu([15 29], [21 24])
%_______________________________________________________________________
% Erling Hugo Jensen, 29/09/09
%
% See also geqPoissonFromVpVs

%% Check input parameters
narginchk(2,2);

%% Calculate 
nu = (3*K - 2*Mu)./(2*(3*K + Mu));
