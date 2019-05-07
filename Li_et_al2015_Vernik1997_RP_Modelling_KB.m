clear; clc; close all;

% Loading measured data froim Vernik 1997
%%% Import the data
[~, ~, raw] = xlsread('M:\matlab\anton\Vernik 1997.xlsx','Appendix A Vernik 1997','A2:J71');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
stringVectors = string(raw(:,1));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[2,3,4,5,6,7,8,9,10]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));
%% Create table
Vernik1 = table;

%% Allocate imported array to column variable names
Vernik1.Formation = categorical(stringVectors(:,1));
Vernik1.Depthft = data(:,1);
Vernik1.Depthm = data(:,2);
Vernik1.Porosity = data(:,3);
Vernik1.BulkDensitygcm3 = data(:,4);    %(Vernik1.Formation=='North Sea')?
Vernik1.Kerogen = data(:,5);
Vernik1.Vp0 = data(:,6);
Vernik1.Vp90 = data(:,7);
Vernik1.Vs0 = data(:,8);
Vernik1.Vsh90 = data(:,9);
Vernik1.AI =(Vernik1.Vp0+Vernik1.Vp90)./2.*Vernik1.BulkDensitygcm3*1000;
Vernik1.SI =(Vernik1.Vs0+Vernik1.Vsh90)./2.*Vernik1.BulkDensitygcm3*1000;
Vernik1.PoissonsRatio0 =(Vernik1.Vp0.^2 - 2*Vernik1.Vs0.^2)./(2*(Vernik1.Vp0.^2-Vernik1.Vs0.^2));
Vernik1.PoissonsRatio90 =(Vernik1.Vp90.^2 - 2*Vernik1.Vsh90.^2)./(2*(Vernik1.Vp90.^2-Vernik1.Vsh90.^2));
Vernik1.PoissonsRatio = (Vernik1.PoissonsRatio0+Vernik1.PoissonsRatio90/2);
Vernik1.Lame = Vernik1.BulkDensitygcm3.*(Vernik1.Vp0.^2-2*Vernik1.Vs0.^2);
Vernik1.Mu = Vernik1.BulkDensitygcm3.*Vernik1.Vs0.^2;
% Defining subplots
n_column_subplots=2;    % columns
n_row_subplots=3;       % rows

FontS = 7 ;
LFS = 15 ;
MarkerS = 6;
Symbol = ['.'];

LW = 1; % Defining line width
color = ['k';'k'; 'k'; 'k'; 'k'; 'k'; 'k']; % ['g';'y'; 'r'; 'c'; 'm'; 'k'; 'y']; % Definerer farverækkefølgen i plots
Symbols = ['^';'o'; 's'; 'x'; 'd'; 'v'; '*']; % Definerer farverækkefølgen i plots
%%
FM = {'Bakken' 'Bazhenov'  'Monterey' 'Niobrara' 'North Sea' 'Japan' 'Lockatong' 'Woodford'};
for i = 1:8

subplot(3,8,i)
Fm = FM(i);
plot(Vernik1.Kerogen(Vernik1.Formation==Fm)*100,...
    Vernik1.Vp0(Vernik1.Formation==Fm),'*')
hold on
plot(Vernik1.Kerogen(Vernik1.Formation==Fm)*100,...
    Vernik1.Vp90(Vernik1.Formation==Fm),'^')
title(Fm)
    xlim([0 max(Vernik1.Kerogen)*100])
    ylim([0 max(Vernik1.Vp90)])

subplot(3,8,i+8)
plot(Vernik1.Kerogen(Vernik1.Formation==Fm)*100,...
    Vernik1.Vs0(Vernik1.Formation==Fm),'*')
hold on
plot(Vernik1.Kerogen(Vernik1.Formation==Fm)*100,...
    Vernik1.Vsh90(Vernik1.Formation==Fm),'^')
title(Fm)
    xlim([0 max(Vernik1.Kerogen)*100])
    ylim([0 max(Vernik1.Vsh90)])

subplot(3,8,i+16)
plot(Vernik1.Kerogen(Vernik1.Formation==Fm)*100,...
    Vernik1.PoissonsRatio0(Vernik1.Formation==Fm),'*')
hold on
plot(Vernik1.Kerogen(Vernik1.Formation==Fm)*100,...
    Vernik1.PoissonsRatio90(Vernik1.Formation==Fm),'^')
title(Fm)
    xlim([0 max(Vernik1.Kerogen)*100])
    ylim([0 max(Vernik1.PoissonsRatio90)])
    
    
    

end
subplot(3,8,1)
ylabel('Vp [km/s]')
legend ('Vp0','Vp90','Location','best' )
subplot(3,8,9)
ylabel('Vs [km/s]')
legend ('Vsh0','Vp90','Location','best' )
subplot(3,8,17)
ylabel('Poissonsratio')
legend ('PR0','PR90','Location','best' )

%% RPM

N = 100;
rock = struct;
rock.Vker = linspace(0,0.3,N);
rock.PHIker = linspace(0,0.5,N);
rock.PHImatrix = 0.05;

% Solid constituents
rock.solid(1).name = 'quartz';
rock.solid(1).vol = 0.5-rock.Vker;
rock.solid(1).k = 38;
rock.solid(1).mu = 44;
rock.solid(1).rho = 2.65;
rock.solid(2).name = 'clay';
rock.solid(2).vol = 0.3;
rock.solid(2).k = 25;
rock.solid(2).mu = 9;
rock.solid(2).rho = 2.5;
rock.solid(3).name = 'kerogen';
rock.solid(3).vol = rock.Vker;
rock.solid(3).k = 6.78;
rock.solid(3).mu = 2.02;
rock.solid(3).rho = 1.4;
rock.solid(4).name = 'dolomite';
rock.solid(4).vol = 0.15;
rock.solid(4).k = 95;
rock.solid(4).mu = 45;
rock.solid(4).rho = 2.87;

% Fluid constituents
rock.fluid(1).name = 'water';
rock.fluid(1).vol = 0;
rock.fluid(1).k = 2.25;
rock.fluid(1).rho = 1.04;
rock.fluid(2).name = 'oil';
rock.fluid(2).vol = 0;
rock.fluid(2).k = 0.57;
rock.fluid(2).rho = 0.7;
rock.fluid(3).name = 'gas';
rock.fluid(3).vol = 1;
rock.fluid(3).k = 0.04;
rock.fluid(3).rho = 0.11;

% Fluid mixture
rock.effFluid.k = geqReuss([rock.fluid(:).k],[rock.fluid(:).vol]);
rock.effFluid.rho = geqEffectiveDensity([rock.fluid(:).rho],[rock.fluid(:).vol]);

% Step 1: Mix kerogen/oil/gas/water using DEM
for n = 1:N
    [rock.ker.kDEM(n), rock.ker.muDEM(n)] = geqDEM([rock.solid(3).k rock.effFluid.k],...
        [rock.solid(3).mu 0],...
        1,... % Aspect ratio of kerogen porosity
        rock.PHIker(n));
    
    rock.ker.rho(n) = (1-rock.PHIker(n))*rock.solid(3).rho + rock.PHIker(n)*rock.effFluid.rho;
    
end

% Step 2: Mix (1) kerogen/oil/gas/water + (2) matrix + (3) shale matrix
% porosity, using SCA

for u = 1:N
    for v = 1:N
        % Berryman's SCA: Bulk and shear modulus
        [rock.kSCA(u,v), rock.muSCA(u,v)] = berryscm(...
            [rock.ker.kDEM(u) rock.solid(1).k rock.solid(2).k rock.solid(4).k rock.fluid(1).k],... % Bulk moduli
            [rock.ker.muDEM(u) rock.solid(1).mu rock.solid(2).mu rock.solid(4).mu 0],... % Shear moduli
            0.1,... % Aspect ratio of shale matrix pore
            [rock.Vker(v) rock.solid(1).vol(v) rock.solid(2).vol rock.solid(4).vol rock.PHImatrix]); % Volue fractions of the components
        
        % Effective density       
        rock.rhoEff(u,v) = rock.ker.rho(u).*rock.Vker(v) + rock.solid(1).rho*rock.solid(1).vol(v) +...
            rock.solid(2).rho*rock.solid(2).vol + rock.solid(4).rho*rock.solid(4).vol +...
            rock.fluid(1).rho*rock.PHImatrix;
    end
end

rock.vpSCA = geqVpFromKMuRho(rock.kSCA,rock.muSCA,rock.rhoEff);
rock.vsSCA = geqVsFromMuRho(rock.muSCA,rock.rhoEff);
rock.aiSCA = rock.vpSCA .* rock.rhoEff;
rock.psSCA = rock.vpSCA ./ rock.vsSCA;
rock.nuSCA = geqPoissonFromKMu(rock.kSCA,rock.muSCA);

% Plot 2D heatmap
figure
subplot(2,2,1)
imagesc('XData',rock.PHIker,'YData',rock.Vker,'CData',rock.vpSCA')
xlabel('\phi_{K} (v/v)'),ylabel('K (v/v)')
set(gca,'XLim',[min(rock.PHIker) max(rock.PHIker)],'YLim',[min(rock.Vker) max(rock.Vker)])
cb1 = colorbar; cb1.Label.String = 'V_p (km/s)';
subplot(2,2,2)
imagesc('XData',rock.PHIker,'YData',rock.Vker,'CData',rock.vsSCA')
xlabel('\phi_{K} (v/v)'),ylabel('K (v/v)')
set(gca,'XLim',[min(rock.PHIker) max(rock.PHIker)],'YLim',[min(rock.Vker) max(rock.Vker)])
cb2 = colorbar; cb2.Label.String = 'V_s (km/s)';
subplot(2,2,3)
imagesc('XData',rock.PHIker,'YData',rock.Vker,'CData',rock.rhoEff')
xlabel('\phi_{K} (v/v)'),ylabel('K (v/v)')
set(gca,'XLim',[min(rock.PHIker) max(rock.PHIker)],'YLim',[min(rock.Vker) max(rock.Vker)])
cb3 = colorbar; cb3.Label.String = '\rho (g/cm^3)';
colormap jet

% Plot curves
figure
subplot(2,2,1)
plot(rock.PHIker,rock.vpSCA(:,1:round((N)^1/3):end))
hold on
scatter(Vernik1.Porosity/1e2,Vernik1.Vp0,20,Vernik1.Kerogen.*100,'filled')
xlabel('\phi_{K} (v/v)'),ylabel('V_p (km/s)')
set(gca,'XLim',[min(rock.PHIker) max(rock.PHIker)])
lg=legend(string(rock.Vker(:,1:round((N)^1/3):end))); lg.Title.String = 'Model: K (v/v)';
subplot(2,2,2)
plot(rock.PHIker,rock.vsSCA(:,1:round((N)^1/3):end))
hold on
scatter(Vernik1.Porosity/1e2,Vernik1.Vs0,20,Vernik1.Kerogen.*100,'filled')
xlabel('\phi_{K} (v/v)'),ylabel('V_s (km/s)')
set(gca,'XLim',[min(rock.PHIker) max(rock.PHIker)])
lg=legend(string(rock.Vker(:,1:round((N)^1/3):end))); lg.Title.String = 'Model: K (v/v)';
cb = colorbar; cb.Label.String = 'Data: K (v/v)';
subplot(2,2,3)
plot(rock.PHIker,rock.nuSCA(:,1:round((N)^1/3):end))
hold on
scatter(Vernik1.Porosity/1e2,Vernik1.PoissonsRatio,20,Vernik1.Kerogen.*100,'filled')
xlabel('\phi_{K} (v/v)'),ylabel('\nu (unitless)')
set(gca,'XLim',[min(rock.PHIker) max(rock.PHIker)])
lg=legend(string(rock.Vker(:,1:round((N)^1/3):end))); lg.Title.String = 'Model: K (v/v)';
colormap(flipud(pmkmp(100,'CubicL')));

