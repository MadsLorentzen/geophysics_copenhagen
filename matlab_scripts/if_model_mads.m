%%% The results match with Figure 9 in <Chalk: composition, diagenesis and physical properties>
clear;

K_a=0.000131;  %% bulk modulus of air

K_ca=71;  
G_ca=30;       %% chalk
Rho_ca=2.71;

K_w=2.2;
G_w=0;         %% water
Rho_w=1.0;


Por=0:0.01:1;
Rho=Rho_ca*(1-Por)+Por*Rho_w;

n=1;
WT = input('wet or dry sample? 0-dry, 1-wet:\n');

for IF=0:0.1:1
	f1=IF*(1-Por);
	f2=Por+(1-IF)*(1-Por);
	
	if WT ==0
	   K_sus=K_a;                                     %%% for dry sample
	elseif WT==1
	   K_sus=(Por./K_w+(1-Por).*(1-IF)/K_ca).^(-1);   %%%  for wet sample
	end

	K=K_ca+f2./((K_sus-K_ca).^(-1)+f1.*(K_ca+4/3*G_ca).^(-1));
	G=G_ca+ f2./(2*f1.*(K_ca+2*G_ca)/(5*G_ca*(K_ca+4/3*G_ca))-1/G_ca);

	M=K+4/3*G;
	Vp=(M/Rho)^0.5;
	Vs=(G/Rho)^0.5;

	M_out(n,:)=M; 
	G_out(n,:)=G; 
	n=n+1;
end

figure; plot(Por,M_out); xlabel('Porosity'); ylabel('M (GPa)'); grid on;
legend('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1');

figure; plot(Por,G_out); xlabel('Porosity'); ylabel('G (GPa)'); grid on;
legend('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1');

%%% The IF=1 bound should be Hashin-Strikman upper-bound


