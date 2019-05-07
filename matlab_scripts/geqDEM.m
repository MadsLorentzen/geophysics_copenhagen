function [k,mu] = geqDEM(K, MU, Aspect, Concentration)
%% function [KEff, MuEff] = geqDEM(K, MU, Aspect, Concentration)
% calculate bulk and shear modulus using DEM
% 
%   KEff            - effective bulk modulus
%   MuEff           - effective shear modulus
%   K               - bulk modulus; [host incl1 incl2 ... inclN]
%   MU              - shear modulus; [host incl1 incl2 ... inclN]
%   Aspect          - aspect ratios; [incl1 incl2 ... inclN]
%   Concentration   - asp.r.conc.; sum([incl1 incl2 ... inclN]) == 1-host
%
% If multiple inclusions are used, then aspect ratio and concentration of
% each must be specified.
% If multiple inclusions with various aspect ratios are used for the same
% inclusion, then moduli of the inclusion must be repeated the same
% number of times as aspect ratios.
% 
% Original code by T. Mukerji, 1997.
% Extended by Erling H. Jensen (2009) to allow for multiple aspect ratios
% and multiple inclusions.
%_______________________________________________________________________
% Examples
% % Porosity = 0.25, spherical pores
% [KEff, MuEff] = geqDEM([39 1], [44 0], [1.0], [0.25])
% % Porosity = 0.25, various shapded pores
% [KEff, MuEff] = geqDEM([39 1 1 1], [44 0 0 0], [1.0 0.3 0.01], [0.20 0.04 0.01])
%_______________________________________________________________________
% Erling Hugo Jensen, 03/12/09

%% Check input parameters
error(nargchk(4,4,nargin,'struct'));

%% The code below is the original, with some minor modifications.

%DEM1 - Effective elastic moduli using Differential Effective Medium
%      formulation. Returns effective moduli at the porosity POR
%      specified in the input.
%
%[K,MU]=DEM1(K1,MU1,K2,MU2,ASP,PHIC,POR)
%
%	K1, MU1:	Bulk and shear moduli of background matrix
%	K2, MU2:	Bulk and shear moduli of inclusions
%	ASP:		Aspect ratio of inclusions
%			<1 for oblate spheroids; >1 for prolate spheroids
%	PHIC:		percolation porosity for modified DEM model
%			=1 for usual DEM
%	K, MU:		Effective bulk and shear moduli
%	POR:		Porosity, fraction of phase 2.
%			For the modified DEM, where phase two is the
%			critical phase, POR is the actual porosity.

%Written by T. Mukerji, 1997

global DEMINPT;
DEMINPT={};

%mu2=3*k2*(1-2*nu2)/(2-2*nu2);

nrOfInclusion           = numel(Aspect);
k                       = 0.0;
mu                      = 0.0;
% tfinal=porosity./phic;

%save deminpt;

%[tout, yout]=ode45m('demyprime',0.00,tfinal,[k1; mu1],1e-10);
% [tout, yout]=ode45m('geqDEMYPrime',0.00,tfinal,[K(1); MU(1)],1e-5);
% for j = 1 : nrOfInclusion
%     DEMINPT(1)=K(1); 
%     DEMINPT(2)=MU(1); 
%     DEMINPT(3)=K(j+1); 
%     DEMINPT(4)=MU(j+1); 
%     DEMINPT(5)=Aspect(j) 
%     DEMINPT(6)=1; 
%     tfinal=Concentration(j)./1;
    DEMINPT.K1              = K(1);
    DEMINPT.MU1             = MU(1);
    DEMINPT.K               = K(2:end); 
    DEMINPT.MU              = MU(2:end); 
    DEMINPT.ASP             = Aspect; 

    % In case of multiple inclusions replacing all of the host, the DEM
    % algorithm gives wrong result. This is avoided by leaving a very small
    % part of the host left, and replacing it with 0.999999 of the
    % inclusion material.
    if (numel(Concentration) > 1) && (sum(Concentration) == 1)
        Concentration = Concentration*0.999999;
    end
    DEMINPT.Concentration   = Concentration;
    
    DEMINPT.Phic            = 1; 
    tfinal=sum(Concentration);
    if (tfinal == 1) 
%         warning('Assumes all inclusions have the same elastic moduli');
        k                       = DEMINPT.K(end);
        mu                      = DEMINPT.MU(end);
    else
%     if (Concentration(j) == 1)
%         k                       = k+DEMINPT(3); 
%         mu                      = mu+DEMINPT(4);  
%     else
%         [tout, yout]             = ...
%            ode45m('geqDEMYPrime',0.00,DEMINPT.Concentration,[K(1); MU(1)],1e-5);
        [tout, yout]             = ...
           ode45m('geqDEMYPrime',0.00,tfinal,[K(1); MU(1)],1e-5);

        n                       = length(tout);
        k                       = k+real(yout(n,1)); 
        mu                      = mu+real(yout(n,2));  
%     end
% end
    end

% kv=real(yout(:,1)); muv=real(yout(:,2));
% porv=phic.*tout;
%if nargout==0
%plot(por,k,'-g',por,mu,'--y', 'linewidth', 1);
%end;
clear DEMINPT;
