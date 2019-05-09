function EffectiveElasticModuli = geqReuss(property1, volFraction1, property2, volFraction2)
%% function EffectiveElasticModuli = geqReuss(ElasticModuliVec, VolumeFractionVec)
% calculates the effective elastic moduli using the Reuss model. (... same
% as Wood's equation for homogeneous mixing of fluids).
% elasticModuli = 1/(sum(VolumeFractionVec ./ ElasticModuliVec));
% 
%   EffectiveElasticModuli       - effective elastic moduli
%   ElasticModuliVec   	- [array]: list of elastic moduli 
%   VolumeFractionVec  	- [array]: corresponding list of volume fractions
% 
% - Use vector notation for single calculation with multiple constituents.
% - Use matrix notation for multiple calculation with multiple constituents
%   where each row is one calculation and the elements on each row
%   correponds to a constituent.
%_______________________________________________________________________
% Examples
% BulkModulus       = geqReuss([37 25], [0.80 0.20]);
% bulkModulus       = geqReuss(37, 0.80, 25, 0.20);     % same as above
% bulkModulus       = geqReuss(37, 0.80, 25);           % same as above
% BulkModulus       = geqReuss([37 25 1.5], [0.80 0.15 0.05]);
% BulkModulus       = geqReuss([37 25 1.5; 38 24 1.4], ...
%                       [0.80 0.15 0.05; 0.8 0.12 0.08])
%_______________________________________________________________________
% Erling Hugo Jensen, 19/11/08
%
% See also geqHill, geqVoigt.

% Revised Erling Hugo Jensen 29/08/11.
%   - mainly extended for handling more than two constituents.


%% Check input parameters
narginchk(2,4);

%% Init variables
ElasticModuliVec            = [];       % Constituent densities
VolumeFractionVec           = [];       % with corresponding volume fractions
EffectiveElasticModuli	    = nan;      % Effective elastic moduli

%% 2 variable input of vectors typically used for more than 2 constituents
isInputArray                = isArray(property1) && isArray(volFraction1);
isInputMatrix               = prod((size(property1)>1)*1.) && prod((size(property1)>1)*1.);
if (nargin == 2) && (isInputArray || isInputMatrix)
    if (sum(volFraction1) > 1)
        warning(['com:Sschwss:Geophysics:' mfilename ':aboveUnitVolumeFraction'], ...
            ['Sum of volume fractions = ' num2str(sum(volFraction1)) ' > 1.0.']);
    end
    
    ElasticModuliVec        = property1;
    VolumeFractionVec       = volFraction1;
    
    %% Calculate effective property
    EffectiveElasticModuli     = 1 ./ (sum(VolumeFractionVec ./ ElasticModuliVec,2));
elseif (nargin > 2)

%% Standard input for two density fractions
    if (nargin < 4)
        volFraction2        = 1 - volFraction1;
    end
    if sum((volFraction1 + volFraction2) > 1) > 0 
        warning(['com:Sschwss:Geophysics:' mfilename ':aboveUnitVolumeFraction'], ...
            ['Sum of volume fractions = ' num2str(volFraction1 + volFraction2) ...
                ' > 1.0. Volume fraction for second constituent has been reduced.']);
        volFraction2     	= 1 - volFraction1;
    end


    %% Calculate effective property
    EffectiveElasticModuli     = 1 ./ (volFraction1 ./ property1 + volFraction2 ./ property2);
else
    
%%  Wrong input
    error(['com:Sschwss:Geophysics:' mfilename ':badInput'], 'Error in input variables');
end

