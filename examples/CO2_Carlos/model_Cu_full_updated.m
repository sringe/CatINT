%Calculation for flat cathode surface concentrations in the electrochemical
%reduction of CO2 in KHCO3 solutions
%------------------------------
%Data acquisition 
%------------------------------
global delta nmesh
fileID1 = fopen('Poly_Cu_JC.txt');
filename = 'Poly_Cu_JC.txt';
metal = fgetl(fileID1);
tline = fgetl(fileID1);
title = strsplit(tline);
npot = dlmread(filename,'', [2 0 2 0]);
npot;
row = npot+5;
nprod = dlmread(filename,'', [3 0 3 0]);
col = nprod+1;
HCO3m_nominal(1,1) = dlmread(filename,'', [4 0 4 0]);
bic = HCO3m_nominal(1,1);
thickness = dlmread(filename,'', [5 0 5 0]);
Data = dlmread(filename,'', [6 0 row col]);
Data
%Data(1:end,1)

fclose(fileID1);
tp = 'Potential';
tj = 'current_density';
t5 = 'H_{2}';
t6 = 'CO';
t7 = 'CH_{4}';
t8 = 'C_{2}H_{4}';
t9 = 'Formate';
t10 = 'Ethanol';
t11 = 'n-Propanol';
t12 = 'Allylalcohol';
t13 = 'Methanol';
t14 = 'Acetate';
t15 = 'Ethyleneglycol';

%------------------------------
%Data classification into appropiate vectors
%------------------------------

for m = 1:nprod+2
    if strcmp(title(m),tp)== 1
        for n = 1:npot
            vpot(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),tj)==1
        for n = 1:npot
            vj(n,1)= Data(n,m);
        end
    end
    if  strcmp(title(m),t5)==1
        for n = 1:npot
            vcefH2(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t6)==1
        for n = 1:npot
            vcefCO(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t7)==1
        for n = 1:npot
            vcefCH4(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t8)==1
        for n = 1:npot
            vcefC2H4(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t9)==1
        for n = 1:npot
            vcefHCOO(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t10)==1
        for n = 1:npot
            vcefetol(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t11)==1
        for n = 1:npot
            vcefpropol(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t12)==1
        for n = 1:npot
            vcefallyl(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t13)==1
        for n = 1:npot
            vcefmetol(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t14)==1
        for n = 1:npot
            vcefacetate(n,1)= Data(n,m);
        end
    end
    if strcmp(title(m),t15)==1
        for n = 1:npot
            vcefethylgly(n,1)= Data(n,m);
        end
    end
end
vpot
% vj;
% vcefH2;
% vcefCO;
% vcefCH4;
% vcefC2H4;
% vcefHCOO;
% vcefetol;
% vcefpropol;
% vcefallyl;
% vcefmetol;
% vcefacetate;
% vcefethylgly;
%------------------------------
%Calculation of total faradaic yield
%------------------------------
for k=1:npot
    var = 0.0;
    for l=3:nprod+2
        vceftot(k,1)= var + Data(k,l);
        var = vceftot(k,1);
    end
    vcefOHm(k,1)=1-vceftot(k,1);
end
% vceftot;
% vcefOHm;

for i=1:npot
%------------------------------
%Deffinition of Constants
%------------------------------
%Experimental constants

cefH2 = vcefH2(i);%current efficiency for hydrogen formation (dimensionless)
cefCO = vcefCO(i); %current efficiency for carbon monoxide form. (dimensionless)
cefHCOO = vcefHCOO(i); %current efficiency for formate formation (dimensionless)
cefCH4 = vcefCH4(i); %current efficiency for methane formation (dimensionless)
cefC2H4 = vcefC2H4(i); %current efficiency for ethylene formation (dimensionless)
cefetol = vcefetol(i); %current efficiency for ethanol formation (dimensionless)
cefpropol = vcefpropol(i); %current efficiency for n-propanol formation (dimensionless)
cefallyl = vcefallyl(i); %current efficiency for allyl alcohol formation (dimensionless)
cefmetol = vcefmetol(i); %current efficiency for methanol formation (dimensionless)
cefacetate = vcefacetate(i); %current efficiency for acetate formation (dimensionless)
cefethylgly = vcefethylgly(i); %current efficiency for ethylene glycol formation (dimensionless)

if vcefOHm(i) <= 0.0
    cefunknown = 0.0; %current efficiency for undetected compounds (prob. glyoxal) (dimensionless)
else
    cefunknown = vcefOHm(i);
end
% cefunknown;
j = -vj(i)*10% current density at the Cu electrode [A/m2]
delta = thickness; % boundary layer thickness [m]
time = 20.0; %total time of experiment [s]
mu = viscosity(bic)*1.0e-003;% viscosity for KHCO3 aqueous solution [Pa.s]
mu0 = 8.90e-004; % viscosity of water at 25C [Pa.s]
global T R F
T = 298.0; %Temperature [K]

P = 1; %Pressure [atm]
%------------------------------
%Physical constants
F = 96485.333; %Faraday's constant [C/mol of e-]
R = 8.314459848; %Gas constant [J/K/mol]
%------------------------------ 
%electrochemical reactions on the electrode
% 2 electron reactions
%   2H2O + 2e- <-> H2 + 2OH-
%   CO2 + H2O + 2e- <-> CO + 2OH-
%   CO2 + H2O + 2e- <-> HCOO- + OH-
% 8 electron reactions
%   CO2 + 6H2O + 8e- <-> CH4 + 8OH-
% 12 electron reactions
%   2CO2 + 8H2O + 12e- <-> C2H4 + 12 OH-
zeffH2 = 2.0; %electrons exchanged in hydrogen formation (dimensionless)
zeffCO = 2.0; %electrons exchanged in carbon monoxide form. (dimensionless)
zeffHCOO = 2.0; %electrons exchanged in formate formation (dimensionless)
zeffCH4 = 8.0; %electrons exchanged in methane formation (dimensionless)
zeffC2H4 = 12.0; %electrons exchanged in ethylene formation (dimensionless)
zeffetol = 12.0; %electrons exchanged in ethanol formation (dimensionless)
zeffallyl = 16.0; %electrons exchanged in allyl alcohol formation (dimensionless)
zeffpropol = 18.0; %electrons exchanged in propanol formation (dimensionless)
zeffmetol = 6.0; %electrons exchanged in methanol formation (dimensionless)
zeffacetate = 8.0; % electrons exchanged in acetate formation (dimensionless)
zeffethylgly = 10.0; % electrons exchanged in ethylene glycol formation (dimensionless)
zeffunknown = 6.0; % assumed electrons exchanged in undetected compound formation (prob. glyoxal) (dimensionless)

%-----------------------------
global K1b k1f k1r K2 k2f k2r
%-----------------------------
% equilibrium reactions  in the CO2-KHCO3 electrolyte
% CO2(aq) + H2O <-> H2CO3 
KH = 2.63e-003; %Equilibrium constant for CO2(aq)-H2CO3 (dimensionless)
% <1% of dissolved CO2 is present as carbonic acid which dissociates very quickly to bicarbonate ion
% H2CO3(aq) <-> HCO3- + H+    K0 = 1.7e-004 M
% thus , the CO2(aq)-HCO3- equilibrium is given by 
% CO2(aq) + H2O <-> HCO3- + H+
K1a = (4.44e-007)*1000.0; % equilibrium constant for CO2(aq)-HCO3- in acid [mol/m3] 
% in basic solutions where pH>7 the CO2(aq)-HCO3- equilibrium is given by 
% CO2(aq) + HO- <-> HCO3- 
K1b = (4.44e007)/1000.0; % equilibrium constant for CO2(aq)-HCO3- in base [1/mol/m3] 
k1f = (5.93e003)/1000.0; %forward rate constant for HCO3-formation [1/mol/(m3 s)]
k1r = k1f/K1b; %reverse rate constant for CO2(aq) formation from HCO3- [1/s]
%The equilibrium constant for the second proton dissociation is given by
% HCO3- <-> (CO3)2- + H+        K2a = 4.66e-011 M
% The bicarbonate ion may be neutralized by OH- generated on the cathode
% HCO3- + OH- <-> (CO3)2- + H2O  
K2 = (4.66e003)/1000.0; % equilibrium constant for HCO3- -(CO3)2- in base [1/mol/m3]
k2f = (1.0e008)/1000.0; % assumed forward rate constant for (CO3)2- formation limited only by OH- diffusion [1/mol/m3/s]
k2r = k2f/K2; % assumed reverse rate constant for HCO3- formation from (CO3)2-[1/s]
% the equilibrium between the bicarbonate and the other CO2 and (CO3)2-
% species in the electrolyte is given by 
% CO2(aq) + (CO3)2- + H2O <-> 2HCO3-
K3 = 9.52e003; % equilibrium constant for (CO3)2- -CO2(aq)-(CO3)2-  (dimensionless)
%------------------------------
%Calculation of initial equilibrium values for CO2, HCO3-, CO32-, OH- and
%pH
%------------------------------
global C_CO2_i C_CO32m_i C_HCO3m_i C_OHm_i C_K_i pH_i
C_CO2_i = 0.03419*P*1000.0; %initial CO2(aq) bulk concentration at t=0 and Pressure P in [mol/m3] units
C_CO32m_i = ((2*bic*1000.0+K3*C_CO2_i)-...
    (sqrt((2*bic*1000.0+K3*C_CO2_i)^2 ...
    -4.0*(bic*1000.0)^2)))/2;  %initial (CO3)2- bulk concentration at t=0 [mol/m3]
% Initial composition of the bulk electrolyte at t=0
C_HCO3m_i = bic*1000.0-C_CO32m_i; %initial HCO3- bulk concentration at t=0 [mol/m3]
C_K_i = bic*1000.0; %initial K+ bulk concentration at t=0 [mol/m3]
C_OHm_i = C_HCO3m_i/K1b/C_CO2_i; %initial OH- bulk concentration at t=0 [mol/m3]
pH_i = 14+log10(C_OHm_i/1000.0); %initial pH


%-----------------------------
%Diffusion coefficients (ref: Cussler, E. L. (1997). Diffusion: Mass Transfer in Fluid Systems (2nd ed.) )
%Diffusion coefficient for HCOO- from Petr Vanýsek: Ionic conductivity and diffusion at infinite dilution, Handbook of Chemistry and Physics
D0_CO2 = 1.91e-009; %Diffusion coefficient of CO2 in water at 25C at infinite dilution [m/s]
D0_CO32m = 9.23e-010; %Diffusion coefficient of (CO3)2- in water at 25C at infinite dilution [m/s]
D0_HCO3m = 1.185e-009; %Diffusion coefficient of HCO3- in water at 25C at infinite dilution [m/s]
D0_OHm = 5.273e-009; %Diffusion coefficient of HCO3- in water at 25C at infinite dilution [m/s]
D0_H2 = 4.50e-009; %Diffusion coefficient of H2 in water at 25C at infinite dilution [m/s]
D0_CO = 2.03e-009; %Diffusion coefficient of CO in water at 25C at infinite dilution [m/s]
D0_CH4 = 1.49e-009; %Diffusion coefficient of CH4 in water at 25C at infinite dilution [m/s]
D0_C2H4 = 1.87e-009; %Diffusion coefficient of C2H4 in water at 25C at infinite dilution [m/s]
D0_HCOOm = 1.454e-009; %Diffusion coefficient of HCOO- in water at 25C at infinite dilution [m/s]
D0_etol = 0.84e-009; %Diffusion coefficient of Ethanol in water at 25C at infinite dilution [m/s]
D0_propol = 1.3e-009; %Diffusion coefficient of propanol in water at 25C at infinite dilution [m/s]
D0_allyl = 1.1e-009; %Diffusion coefficient of allyl alcohol in water at 25C at infinite dilution [m/s]
D0_metol = 0.84e-009; %Diffusion coefficient of methanol in water at 25C at infinite dilution [m/s]
D0_acetate = 1.089e-009; %Diffusion coefficient of acetate in water at 25C at infinite dilution [m/s]
D0_ethylgly = 1.102e-009; %Diffusion coefficient of ethyleneglycol in water at 25C at 0.02 fraction of dlycol[m/s]
D0_K = 1.957e-009; %Diffusion coefficient of potasium ions in water at 25C at infinite dilution [m/s]
%-----------------------------
%Diffusion coefficients corrected for varying electrolyte concentration
%using Stokes-Einsteinn's equation
global D_CO2 D_CO32m D_HCO3m D_OHm D_H2 D_CO D_CH4 D_C2H4 D_HCOOm D_etol
global D_propol D_allyl D_metol D_acetate D_ethylgly D_K
D_CO2 = D0_CO2*mu0/mu %Diffusion coefficient of CO2 in HCO3- solution at 25C [m/s]
D_CO32m = D0_CO32m*mu0/mu %Diffusion coefficient of (CO3)2- in HCO3- solution at 25C [m/s]
D_HCO3m = D0_HCO3m*mu0/mu; %Diffusion coefficient of HCO3- in HCO3- solution at 25C [m/s]
D_OHm = D0_OHm*mu0/mu; %Diffusion coefficient of OH- in HCO3- solution at 25C [m/s]
D_H2 = D0_H2*mu0/mu; %Diffusion coefficient of H2 in HCO3- solution at 25C [m/s]
D_CO = D0_CO*mu0/mu; %Diffusion coefficient of CO in HCO3- solution at 25C [m/s]
D_CH4 = D0_CH4*mu0/mu; %Diffusion coefficient of CH4 in HCO3- solution at 25C [m/s]
D_C2H4 = D0_C2H4*mu0/mu; %Diffusion coefficient of C2H4 in HCO3- solution at 25C [m/s]
D_HCOOm = D0_HCOOm*mu0/mu; %Diffusion coefficient of C2H4 in HCO3- solution at 25C [m/s]
D_etol = D0_etol*mu0/mu; %Diffusion coefficient of ethanol in HCO3- solution at 25C [m/s]
D_propol = D0_propol*mu0/mu; %Diffusion coefficient of propanol in HCO3- solution at 25C [m/s]
D_allyl = D0_allyl*mu0/mu; %Diffusion coefficient of allyl alcohol in HCO3- solution at 25C [m/s]
D_metol = D0_metol*mu0/mu; %Diffusion coefficient of methanol in HCO3- solution at 25C [m/s]
D_acetate = D0_acetate*mu0/mu; %Diffusion coefficient of acetate in HCO3- solution at 25C [m/s]
D_ethylgly = D0_ethylgly*mu0/mu; %Diffusion coefficient of ethyleneglycol in HCO3- solution at 25C [m/s]
D_K = D0_K*mu0/mu; %Diffusion coefficient of ethyleneglycol in HCO3- solution at 25C [m/s]

% %Ionic mobilities corrected for varying electrolyte concentration
% u_CO32m = D_CO32m %mobility of (CO3)2- in HCO3- solution at 25C [m/s]
% u_HCO3m = D_HCO3m; %mobility of HCO3- in HCO3- solution at 25C [m/s]
% u_OHm = D_OHm; %mobility of OH- in HCO3- solution at 25C [m/s]
% u_HCOOm = D_HCOOm; %mobility of formate in HCO3- solution at 25C [m/s]
% u_acetate = D_acetate; %mobility of acetate in HCO3- solution at 25C [m/s]
% u_K = D_K %mobility of potassium ion in HCO3- solution at 25C [m/s]
%-----------------------------
global CO2consumption OHmformation H2formation COformation CH4formation vpot_i
global C2H4formation HCOOmformation etolformation propolformation 
global allylformation metolformation acetateformation ethylglyformation unknownformation
% Calculation of CO2 consumption rate at cathode surface (mol/m2/s)
CO2consumption =((j/F)*(cefHCOO/zeffHCOO + cefCO/zeffCO + cefCH4/zeffCH4 +...
    + 2.0*cefC2H4/zeffC2H4 + 2.0*(cefetol/zeffetol)+ 3.0*(cefpropol/zeffpropol)...
    + 3.0*(cefallyl/zeffallyl)+ cefmetol/zeffmetol + 2.0*(cefacetate/zeffacetate)...
    + 2.0*(cefethylgly/zeffethylgly)+ 2.0*(cefunknown/zeffunknown) ));
% Calculation of OH- rate of formation at cathode surface (mol/m2/s)
OHmformation = ((j/F)*(cefHCOO/zeffHCOO + 2.0*(cefCO/zeffCO) + ...
    2.0*(cefH2/zeffH2) + 8.0*(cefCH4/zeffCH4) + 12.0*(cefC2H4/zeffC2H4) ...
    + 12.0*(cefetol/zeffetol)+ 16.0*(cefallyl/zeffallyl)+ ...
    18.0*(cefpropol/zeffpropol)+6.0*(cefmetol/zeffmetol)+ ...
    8.0*(cefacetate/zeffacetate)+10.0*(cefethylgly/zeffethylgly)+ ...
    6.0*(cefunknown/zeffunknown)));
% Calculation of H2 rate of formation at cathode surface (mol/m2/s)
H2formation = (j/F)*(cefH2/zeffH2);
% Calculation of CO rate of formation at cathode surface (mol/m2/s)
COformation = (j/F)*(cefCO/zeffCO);
% Calculation of CO rate of formation at cathode surface (mol/m2/s)
CH4formation = (j/F)*(cefCH4/zeffCH4);
% Calculation of CO rate of formation at cathode surface (mol/m2/s)
C2H4formation = (j/F)*(cefC2H4/zeffC2H4);
% Calculation of HCOO- rate of formation at cathode surface (mol/m2/s)
HCOOmformation = (j/F)*(cefHCOO/zeffHCOO);
% Calculation of ethanol rate of formation at cathode surface (mol/m2/s)
etolformation = (j/F)*(cefetol/zeffetol);
% Calculation of propanol rate of formation at cathode surface (mol/m2/s)
propolformation = (j/F)*(cefpropol/zeffpropol);
% Calculation of allyl alcohol rate of formation at cathode surface (mol/m2/s)
allylformation = (j/F)*(cefallyl/zeffallyl);
% Calculation of methanol rate of formation at cathode surface (mol/m2/s)
metolformation = (j/F)*(cefmetol/zeffmetol);
% Calculation of unknown rate of formation at cathode surface (mol/m2/s)
acetateformation = (j/F)*(cefacetate/zeffacetate);
% Calculation of unknown rate of formation at cathode surface (mol/m2/s)
ethylglyformation = (j/F)*(cefethylgly/zeffethylgly);
% Calculation of unknown rate of formation at cathode surface (mol/m2/s)
unknownformation = (j/F)*(cefunknown/zeffunknown);
vpot_i = vpot(i);



% %-----------------------------
% %Diffusion coefficients file
% %-----------------------------
% fileID2 = fopen('temporal_Dif.txt','w');
% fprintf(fileID2,'%e\n',D_CO2);
% fprintf(fileID2,'%e\n',D_CO32m);
% fprintf(fileID2,'%e\n',D_HCO3m);
% fprintf(fileID2,'%e\n',D_OHm);
% fprintf(fileID2,'%e\n',D_H2);
% fprintf(fileID2,'%e\n',D_CO);
% fprintf(fileID2,'%e\n',D_CH4);
% fprintf(fileID2,'%e\n',D_C2H4);
% fprintf(fileID2,'%e\n',D_HCOOm);
% fprintf(fileID2,'%e\n',D_etol);
% fprintf(fileID2,'%e\n',D_propol);
% fprintf(fileID2,'%e\n',D_allyl);
% fprintf(fileID2,'%e\n',D_metol);
% fclose(fileID2);

%-----------------------------
%SOLUTION BLOCK
%-----------------------------

nmesh = 101; % intial mesh
xplot = nmesh; % meshes for plotting the result on x.
tplot = nmesh; % meshes for plotting the result on t.
% solution block.
m = 0; 
x = linspace(0,delta,nmesh); 
t = linspace(0,time,nmesh);

sol = pdepe(m,@pdex5pde,@pdex5ic,@pdex5bc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);
u5 = sol(:,:,5);
u6 = sol(:,:,6);
u7 = sol(:,:,7);
u8 = sol(:,:,8);
u9 = sol(:,:,9);
u10 = sol(:,:,10);
u11 = sol(:,:,11);
u12 = sol(:,:,12);
u13 = sol(:,:,13);
u14 = sol(:,:,14);
u15 = sol(:,:,15);
u16 = sol(:,:,16);
u17 = sol(:,:,17);
%0.05916*(14+log10(abs(sol(:,:,4)/1000))-pH_i)

CO2conc = sol(nmesh,:,1)/1000.0;
HCO3mconc = sol(nmesh,:,2)/1000.0;
CO32mconc = sol(nmesh,:,3)/1000.0;
OHmconc = sol(nmesh,:,4)/1000.0;
pH = 14+log10(abs(OHmconc));

H2conc = sol(nmesh,:,5)/1000.0;
COconc = sol(nmesh,:,6)/1000.0;
CH4conc = sol(nmesh,:,7)/1000.0;
C2H4conc = sol(nmesh,:,8)/1000.0;
HCOOmconc = sol(nmesh,:,9)/1000.0;
etolconc = sol(nmesh,:,10)/1000.0;
propolconc = sol(nmesh,:,11)/1000.0;
allylconc = sol(nmesh,:,12)/1000.0;
metolconc = sol(nmesh,:,13)/1000.0;
acetateconc = sol(nmesh,:,14)/1000.0;
ethylglyconc = sol(nmesh,:,15)/1000.0;
Kconc = sol(nmesh,:,16)/1000.0;
charge_total= (u16-u2-2*u3-u4-u9-u14)/1000.0;
final_charge =Kconc-HCO3mconc-2*CO32mconc-OHmconc-HCOOmconc-acetateconc;
Potential = vpot_i+0.05916*(14+log10(abs(u4)/1000.0)-pH_i);
final_potential = vpot_i+0.05916*(14+log10(OHmconc)-pH_i);
% DeltaP = sol(nmesh,:,17);


% figure;
% surf(x,t,charge_total,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'zlim',[-0.1 0.1]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('charge neutrality[mol/L]');
% view(30,20);
% 
% figure;
% surf(x,t,u17,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('Potential[]');
% view(30,20);
% 
% figure;
% surf(x,t,u16/1000,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('K+[M]');
% view(30,20);
% 
% for l=1:nmesh-1
% deltatime(l)=t(l+1)-t(l);
% deltax(l)=x(l+1)-x(l);
% Derivative1(l,:)=(u17(l+1,:)-u17(l,:))/deltatime(l);
% Derivativex1(:,l)=(u17(:,l+1)-u17(:,l))/deltax(l);
% end
% Derivative1(nmesh,:)=Derivative1(nmesh-1,:);
% Derivativex1(:,nmesh)=Derivativex1(:,nmesh-1);
% deltatime(nmesh)=deltatime(nmesh-1);
% deltax(nmesh)=deltax(nmesh-1);
% epsi = 8.854187817E-12;
% for l=1:nmesh-1
% Derivative2(l,:)=D_K*F/R/T*(Derivative1(l+1,:)-Derivative1(l,:))/...
%     deltatime(l)/epsi;
% Derivativex2(:,l)=(Derivativex1(:,l+1)-Derivative1(:,l))/deltax(l);
% end
% Derivative2(nmesh,:)=Derivative2(nmesh-1,:);
% Derivativex2(:,nmesh)=Derivativex2(:,nmesh-1);


% Derivative2(:,:,1)=u17(:,:,1);
% for i=1:3
%         Derivative2(:,:,i+1)=u17(:,:,i)
% % for i=1:nmesh-1
% %     Derivative = Derivative1-Derivative2(i)
% end
% figure;
% surf(x,t,Derivative1,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('First time derivative');
% view(30,20);
% 
% figure;
% surf(x,t,Derivative2,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('second time derivative');
% view(30,20);
% 
% figure;
% surf(x,t,Derivativex1,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('First x derivative');
% view(30,20);
% 
% figure;
% surf(x,t,Derivativex2,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('second x derivative');
% view(30,20);

% figure;
% surf(x,t,u17,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('Electric Potential [V]');
% view(30,20);

% figure;
% surf(x,t,Potential,'edgecolor','none');
% set(gca,'xlim',[0.0 0.8e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('Nernstian Potential [V]');
% view(30,20);

% figure;
% surf(x,t,u1/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('CO_{2} concentration [M]');
% view(30,20);
% 
% figure;
% surf(x,t,u2/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('HCO_{3}^{-} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('HCO_{3}^{-} concentration [M]');
% view(30,20);
% 
% figure;
% surf(x,t,u3/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('CO_{3}^{2-} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('CO_{3}^{2-} concentration [M]');
% 
% figure;
% surf(x,t,14+log10(abs(u4/1000.0)),'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('pH (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('pH');
% 
% figure;
% surf(x,t,u5/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('H_{2} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('H_{2} concentration [M]');
% 
% figure;
% surf(x,t,u6/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('CO (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('CO concentration [M]');
% 
% figure;
% surf(x,t,u7/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('CH_{4} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('CH_{4} concentration [M]');
% 
% figure;
% surf(x,t,u8/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('C_{2}H_{4} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('C_{2}H_{4} concentration [M]');
% 
% figure;
% surf(x,t,u9/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('HCOO^{-} (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('HCOO^{-} concentration [M]');
% 
% figure;
% surf(x,t,u10/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('Ethanol (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('Ethanol concentration [M]');
% 
% figure;
% surf(x,t,u11/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('n-Propanol (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('n-Propanol concentration [M]');
% 
% figure;
% surf(x,t,u12/1000.0,'edgecolor','none');
% set(gca,'xlim',[0.0 1e-004]);
% set(gca,'fontsize',11.5);
% %title('Allyl alcohol (x,t) in KHCO_{3} = 0.1 M, \delta = 0.01 cm');
% xlabel('Distance x [m]');
% ylabel('Time t [s]');
% zlabel('Allyl alcohol concentration [M]');

figure
%title({'Steady state CO_{2}, HCOO^{-}, OH^{-} conc. in KHCO_{3} = 0.1 M,\delta = 0.01 cm';' '});
hold on

yyaxis right
plot(x,OHmconc,'LineWidth',1.5);
ylabel('OH^{-} concentration [M]');
set(gca,'xlim',[0.0 delta]);
set(gca,'fontsize',11.5);
xlabel('Distance x [m]');

yyaxis left
plot(x,CO2conc,x,HCO3mconc,x,CO32mconc,x,Kconc,'LineWidth',1.5);
ylabel('Concentration [M]');

hold off
legend('CO_{2}','HCO_{3}^{-}','CO_{3}^{2-}','K^{+}','OH^{-}','Location','west');

figure
plot(x,H2conc,x,COconc,x,CH4conc,x,C2H4conc,x,HCOOmconc,x,etolconc,...
    x,propolconc,x,allylconc,x,metolconc,'LineWidth',1.5);
set(gca,'xlim',[0.0 delta]);
set(gca,'fontsize',11.5);
xlabel('Distance x [m]');
ylabel('Product concentration [M]');
legend('H_{2}','CO','CH_{4}','C_{2}H_{4}','HCOO^{-}','Ethanol','n-Propanol','Allyl alcohol','Methanol','Location','northwest');
%title({'Products conc. in KHCO_{3} = 0.1 M, \delta = 0.01 cm';' '});

    figure
    plot(x,pH,'LineWidth',1.5);
    set(gca,'xlim',[0.0 delta]);
    set(gca,'fontsize',11.5);
    %set(gca,'linewidth',1);
    xlabel('Distance x [m]');
    ylabel('pH');
    legend('5 mA cm^{-2}','Location','northwest');
    %title({'Steady state pH in KHCO_{3} = 0.1 M, \delta = 0.01 cm';' '});


%final concentrations on the electrode surface
u1final = u1(nmesh,nmesh)/1000.0;
u2final = u2(nmesh,nmesh)/1000.0;
u3final = u3(nmesh,nmesh)/1000.0;
u4final = u4(nmesh,nmesh)/1000.0;
u5final = u5(nmesh,nmesh)/1000.0;
u6final = u6(nmesh,nmesh)/1000.0;
u7final = u7(nmesh,nmesh)/1000.0;
u8final = u8(nmesh,nmesh)/1000.0;
u9final = u9(nmesh,nmesh)/1000.0;
u10final = u10(nmesh,nmesh)/1000.0;
u11final = u11(nmesh,nmesh)/1000.0;
u12final = u12(nmesh,nmesh)/1000.0;
u13final = u13(nmesh,nmesh)/1000.0;
u14final = u14(nmesh,nmesh)/1000.0;
u15final = u15(nmesh,nmesh)/1000.0;
u16final = u16(nmesh,nmesh)/1000.0;
u17final = u17(nmesh,nmesh);
pHfinal = 14 + log10(u4final);
str = string({'results_multipot_Cu_full_model',i,'.txt'});
newstr = join(str,'_');
fileID4 = fopen(newstr,'w');
fprintf(fileID4,'%8s %e\n','Potent.',vpot(i));
fprintf(fileID4,'%8s %e\n','Current',vj(i));
fprintf(fileID4,'%8s %12s %12s\n','Product','Concen.[M]','Current_eff');
fprintf(fileID4,'%8s %e %e\n','CO2', u1final, 0);
fprintf(fileID4,'%8s %e %e\n','HCO3m', u2final, 0);
fprintf(fileID4,'%8s %e %e\n','CO32m', u3final, 0);
fprintf(fileID4,'%8s %e %e\n','OHm', u4final, 0);
fprintf(fileID4,'%8s %e %e\n','H2', u5final, cefH2);
fprintf(fileID4,'%8s %e %e\n','CO', u6final, cefCO);
fprintf(fileID4,'%8s %e %e\n','CH4', u7final, cefCH4);
fprintf(fileID4,'%8s %e %e\n','C2H4', u8final, cefC2H4);
fprintf(fileID4,'%8s %e %e\n','HCOOm', u9final, cefHCOO);
fprintf(fileID4,'%8s %e %e\n','Ethol', u10final, cefetol);
fprintf(fileID4,'%8s %e %e\n','Propol', u11final, cefpropol);
fprintf(fileID4,'%8s %e %e\n','Allylol', u12final, cefallyl);
fprintf(fileID4,'%8s %e %e\n','Metol', u13final, cefmetol);
fprintf(fileID4,'%8s %e %e\n','Acetate', u14final,cefacetate);
fprintf(fileID4,'%8s %e %e\n','Ethylgly', u15final,cefethylgly);
fprintf(fileID4,'%8s %e %e\n','K', u16final,0);
fprintf(fileID4,'%8s %e %e\n','Electric_potential_surface', u17final,0);
fclose(fileID4);
end

function mu = viscosity(x)
Data = dlmread('viscosity.txt','', [0 0 18 1]);
if x<Data(2,1) | x>Data(19,1)
    mu = 'viscosity outside limits';
else
    for n=1:18
    if x>=Data(n,1) & x<=Data(n+1,1)
       mu = Data(n,2)+(x-Data(n,1))*((Data(n+1,2)-Data(n,2))/(Data(n+1,1)-Data(n,1)));
    end
    end
end
end

function [c,f,s] = pdex5pde(x,t,u,DuDx)
    %Diffusion coefficients corrected for varying electrolyte concentration
    %using Stokes-Einsteinn's equation

	global D_CO2 D_CO32m D_HCO3m D_OHm D_H2 D_CO D_CH4 D_C2H4 D_HCOOm D_etol
	global D_propol D_allyl D_metol D_acetate D_ethylgly D_K delta nmesh tstart
	global K1b k1f k1r K2 k2f k2r pH_i vpot_i OHmstart tendb T R F e epsi avog 
	global f1
	e = 1.602176620898E-19;
	epsi = 8.854187817E-12;
	avog = 6.02214085E23;
	% syms y(x)
	% [V] = odeToVectorField(diff(y, 2) == -e*(u(2)+2.0*u(3)+u(4)+u(9)+u(14)...
	%     -u(16))*avog/epsi/(78.304));
	% M = matlabFunction(V,'vars', {'x','Y'})
	% sol = ode45(M,x,Y);
	%-----------------------------
	c = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1E-80];
	f = [D_CO2*DuDx(1); D_HCO3m*DuDx(2)-D_HCO3m*1.0*u(2)*F/R/T*DuDx(17);...
	    D_CO32m*DuDx(3)-D_CO32m*4.0*u(3)*F/R/T*DuDx(17); ...
	    D_OHm*DuDx(4)-D_OHm*1.0*u(4)*F/R/T*DuDx(17); ...
	    D_H2*DuDx(5); D_CO*DuDx(6); D_CH4*DuDx(7); D_C2H4*DuDx(8); ...
	    D_HCOOm*(DuDx(9)-1.0*u(9)*F/R/T*DuDx(17));D_etol*DuDx(10);...
	    D_propol*DuDx(11); D_allyl*DuDx(12);...
	    D_metol*DuDx(13); D_acetate*(DuDx(14)-1.0*u(14)*F/R/T*DuDx(17));...
	    D_ethylgly*DuDx(15);...
	    D_K*DuDx(16)+D_K*1.0*u(16)*F/R/T*DuDx(17); DuDx(17)];
	s1 = -k1f*u(1)*u(4) + k1r*u(2);
	s2 = k1f*u(1)*u(4) - k1r*u(2) - k2f*u(2)*u(4) + k2r*u(3);
	s3 = k2f*u(2)*u(4)-k2r*u(3);
	s4 = -k1f*u(1)*u(4) + k1r*u(2) - k2f*u(2)*u(4) + k2r*u(3);
	s5 = 0;
	s6 = 0;
	s7 = 0;
	s8 = 0;
	s9 = 0;
	s10 = 0;
	s11 = 0;
s12 = 0;
s13 = 0;
s14 = 0;
s15 = 0;
s16 = 0;
s17 = -e*(u(2)+2.0*u(3)+u(4)+u(9)+u(14)-u(16))*avog/epsi/(78.304);
s = [s1; s2; s3; s4; s5; s6; s7; s8; s9; s10; s11; s12; s13; s14; s15;...
    s16; s17];
end
%u(16)-u(2)-2.0*u(3)-u(4)-u(9)-u(14)
% --------------------------------------------------------------------------

function u0 = pdex5ic(x)
% Initial composition of the bulk electrolyte at t=0
global C_CO2_i C_CO32m_i C_HCO3m_i C_OHm_i C_K_i pH_i vpot_i
e = 1.602176620898E-19;
C_H2_i = 0.0*1000.0; %initial H2 bulk concentration at t=0 [mol/m3]
C_CO_i = 0.0*1000.0; %initial CO bulk concentration at t=0 [mol/m3]
C_CH4_i = 0.0*1000.0; %initial CH4 bulk concentration at t=0 [mol/m3]
C_C2H4_i = 0.0*1000.0; %initial C2H4 bulk concentration at t=0 [mol/m3]
C_HCOOm_i = 0.0*1000.0; %initial HCOO- bulk concentration at t=0 [mol/m3]
C_etol_i = 0.0*1000.0; %initial ethanol bulk concentration at t=0 [mol/m3]
C_propol_i = 0.0*1000.0; %initial n-propanol bulk concentration at t=0 [mol/m3]
C_allyl_i = 0.0*1000.0; %initial allyl alcohol bulk concentration at t=0 [mol/m3]
C_metol_i = 0.0*1000.0; %initial methanol bulk concentration at t=0 [mol/m3]
C_acetate_i = 0.0*1000.0; %initial acetate bulk concentration at t=0 [mol/m3]
C_ethylgly_i = 0.0*1000.0; %initial ethylene glycol bulk concentration at t=0 [mol/m3]
%------------------------------
u0 = [C_CO2_i; C_HCO3m_i; C_CO32m_i; C_OHm_i; C_H2_i; C_CO_i; ...
    C_CH4_i; C_C2H4_i; C_HCOOm_i; C_etol_i; C_propol_i; C_allyl_i; ...
    C_metol_i; C_acetate_i; C_ethylgly_i; C_K_i; 0.0 ];
end

% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = pdex5bc(xl,ul,xr,ur,t)
% Initial composition of the bulk electrolyte at t=0
global C_CO2_i C_CO32m_i C_HCO3m_i C_OHm_i C_K_i pH_i vpot_i
e = 1.602176620898E-19;
C_H2_i = 0.0*1000.0; %initial H2 bulk concentration at t=0 [mol/m3]
C_CO_i = 0.0*1000.0; %initial CO bulk concentration at t=0 [mol/m3]
C_CH4_i = 0.0*1000.0; %initial CH4 bulk concentration at t=0 [mol/m3]
C_C2H4_i = 0.0*1000.0; %initial CH4 bulk concentration at t=0 [mol/m3]
C_HCOOm_i = 0.0*1000.0; %initial HCOO- bulk concentration at t=0 [mol/m3]
C_etol_i = 0.0*1000.0; %initial HCOO- bulk concentration at t=0 [mol/m3]
C_propol_i = 0.0*1000.0; %initial n-propanol bulk concentration at t=0 [mol/m3]
C_allyl_i = 0.0*1000.0; %initial allyl alcohol bulk concentration at t=0 [mol/m3]
C_metol_i = 0.0*1000.0; %initial methanol bulk concentration at t=0 [mol/m3]
C_acetate_i = 0.0*1000.0; %initial acetate bulk concentration at t=0 [mol/m3]
C_ethylgly_i = 0.0*1000.0; %initial ethylene glycol bulk concentration at t=0 [mol/m3]
global CO2consumption OHmformation H2formation COformation CH4formation
global C2H4formation HCOOmformation etolformation propolformation 
global allylformation metolformation acetateformation ethylglyformation

%------------------------------
pl = [ul(1)-C_CO2_i; ul(2)-C_HCO3m_i; ul(3)-C_CO32m_i; ul(4)-C_OHm_i; ...
    ul(5)-C_H2_i; ul(6)-C_CO_i; ul(7)-C_CH4_i; ul(8)-C_C2H4_i; ul(9)-C_HCOOm_i;...
    ul(10)-C_etol_i; ul(11)-C_propol_i; ul(12)-C_allyl_i; ul(13)-C_metol_i;...
    ul(14)-C_acetate_i; ul(15)-C_ethylgly_i; ul(16)-C_K_i; ul(17)];
ql = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% if(t<100)
% pr = [CO2consumption*t/100; 0.0; 0.0;...
%     -OHmformation*t/100;H2formation*t/100; -COformation*t/100;...
%     -CH4formation*t/100;-C2H4formation*t/100;...
%      -HCOOmformation*t/100; -etolformation*t/100;...
%     -propolformation*t/100; -allylformation*t/100; -metolformation*t/100;...
%     -acetateformation*t/100;-ethylglyformation*t/100;0.0; ur(17)-vpot_i-...
%     0.05916*(14+log10(abs(ur(4))/1000.0)-pH_i);0.0];
% else
pr = [CO2consumption; 2*H2formation; 0.0; -OHmformation; -H2formation; -COformation;...
    -CH4formation; -C2H4formation; -HCOOmformation; -etolformation;...
    -propolformation; -allylformation; -metolformation; -acetateformation;...
    -ethylglyformation; 0.0; 0.0];

% end
qr = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
end
%-vpot_i-0.05916*(14+log10(ul(4)/1000)-pH_i)
% ur(17)-vpot_i-0.05916*(14+log10(abs(ur(4))/1000.0)-pH_i)