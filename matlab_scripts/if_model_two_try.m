


%%% The results match well with Figure 9 in <Chalk: composition, diagenesis and physical properties>
clear;

K_a=0.000131;
Por_critical=0.7;         %%% critical porosity of chalk(handbook) ---would be 70% based on Ida's paper

K_ca=71;  G_ca=30;         %%% from Ida's paper
%K_ca=70.2; G_ca=29.0;      %%% from rock physics handbook
Rho_ca=2.71;

K_w=2.2;
G_w=0;
Rho_w=1.0;


Por=0:0.01:1;
Rho=Rho_ca*(1-Por)+Por*Rho_w;
n=1;

for IF=0:0.1:1
	f1=IF*(1-Por);
	f2=Por+(1-IF)*(1-Por);

	% f1=IF*(Por_critical-Por)/Por_critical;  %%% using critical porosity as endpoint
	% f2=(Por+(1-IF)*(Por_critical-Por))/Por_critical;
	
	K_sus=K_a;                  %%% true only for dry sample
	K_sus=(Por./K_w+(1-Por).*(1-IF)/K_ca).^(-1);   %%%  for wet sample


	% K=K_ca+f2./((K_sus-K_ca).^(-1)+f1.*(K_ca+4/3*G_ca).^(-1));  %% from Casper Olsen thesis/ Ida 2003/ Tuhin Bhakta
	% G=G_ca+ f2./(2*f1.*(K_ca+2*G_ca)/(5*G_ca*(K_ca+4/3*G_ca))-1/G_ca);

	
	%zeta=G_ca/6*((9*K_ca+8*G_ca)/(K_ca+2*G_ca));
	% K=((Por+(1-IF).*(1-Por))./(K_sus+4/3*G_ca)+(IF.*(1-Por)./(K_ca+4/3*G_ca))).^(-1); %% original Ida in <Chalk: composition, diagenesis and physical properties>
	% G=((Por+(1-IF).*(1-Por))./zeta+IF.*(1-Por)./(G_ca+zeta)).^(-1);                   %% original Ida
	
	zeta=G_ca/6*((9*K_ca+8*G_ca)/(K_ca+2*G_ca));	
	K=((Por+(1-IF).*(1-Por))./(K_sus+4/3*G_ca)+(IF.*(1-Por)./(K_ca+4/3*G_ca))).^(-1)-40; %% calibrated Ida -perfect, 40 replaced as 4/3*G_ca?
	G=((Por+(1-IF).*(1-Por))./zeta+IF.*(1-Por)./(G_ca+zeta)).^(-1)-zeta;                 %% calibrated Ida -perfect, zeta=33.55
	
	
	M=K+4/3*G;
	Vp=(M/Rho)^0.5;
	Vs=(G/Rho)^0.5;

	M_out(n,:)=M; G_out(n,:)=G; n=n+1;
end

figure; plot(Por,M_out); xlabel('Porosity'); ylabel('M (GPa)'); grid on;
legend('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1');

figure; plot(Por,G_out); xlabel('Porosity'); ylabel('G (GPa)'); grid on;
%%% The IF=1 bound should be Hashin-Strikman upper-bound


