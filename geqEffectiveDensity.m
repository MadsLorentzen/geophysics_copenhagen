function density = geqEffectiveDensity(rho1, fraction1, rho2, fraction2)
%% function density = geqEffectiveDensity(DensityVec, VolumeFractionVec)
% calculates the effective density (using arithmetic average). 
% density = sum(DensityVec .* VolumeFractionVec);
%
%   density             - effective density
%   DensityVec          - [array]: list of densities 
%   VolumeFractionVec  	- [array]: corresponding list of volume fractions
% 
%_______________________________________________________________________
% Examples
% density   = geqEffectiveDensity([2.65 2.6], [0.6 0.4]);
% density   = geqEffectiveDensity(2.65, 0.6, 2.6, 0.4);     % same as above
% density   = geqEffectiveDensity(2.65, 0.6, 2.6);          % same as above
% density   = geqEffectiveDensity([2.65 2.6 1.0], [0.5 0.4 0.1]);
%_______________________________________________________________________
% Erling Hugo Jensen, 24/11/08

% Revised Erling Hugo Jensen 29/08/11.
%   - mainly extended for handling more than two constituents.
%
% Revised Erling Hugo Jensen 16/08/13.
%   - modified to handle simultaneous multiple calculations for more than
%   two constituents

%% Check input parameters
narginchk(2,4);

%% Init variables
DensityVec                  = [];       % Constituent densities
VolumeFractionVec           = [];       % with corresponding volume fractions
density                     = nan;      % Effective density

%% 2 variable input of vectors typically used for more than 2 constituents
if (nargin == 2) && (~isAtom(rho1)) && (~isAtom(fraction1))
    [DensityVec VolumeFractionVec]    = rectify(rho1, fraction1);
    if (sum(sum(VolumeFractionVec,2) > 1) > 0)
        warning(['com:Sschwss:Geophysics:' mfilename ':aboveUnitVolumeFraction'], ...
                ['Sum of volume fractions = ' num2str(sum(fraction1)) ' > 1.0.']);
    end
       
    %% Calculate effective property
    density     = sum(DensityVec .* VolumeFractionVec,2);
elseif (nargin > 2)

%% Standard input for two density fractions
    if (nargin < 4)
        fraction2           = 1 - fraction1;
    end
    if sum((fraction1 + fraction2) > 1) > 0     
        warning(['com:Sschwss:Geophysics:' mfilename ':aboveUnitVolumeFraction'], ...
            ['Sum of volume fractions = ' num2str(fraction1 + fraction2) ...
                ' > 1.0. Volume fraction for second constituent has been reduced.']);
        fraction2           = 1 - fraction1;
    end

    DensityVec              = [rho1 rho2];
    VolumeFractionVec       = [fraction1 fraction2];
    
    %% Calculate effective property
    density     = rho1.*fraction1 + rho2.*fraction2;
else
    
%%  Wrong input
    error(['com:Sschwss:Geophysics:' mfilename ':badInput'], 'Error in input variables');
end

