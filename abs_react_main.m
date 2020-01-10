% ABS_REACT
% CH4+CO2 (g), MDEA+H2O (l), CO2 transfer with reaction in liquid bulk
clear all; clc;
T       = 40+273.15;     % K
p       = 7.1*10^(6);    % Pa
wMDEA   = 0.5;           % weight fract of MDEA in liq
u       = 0.1436;        % superficial velo. m/s
S       = 0.1;           % submergence of the liquid m
nLin    = 52.82493556;   % inlet liquid flow mol/s
nVin    = 347.0251241;   % inlet gas flow mol/s
xCO2in  = 0;             % inlet conc of CO2 in liq
yCO2in  = 0.05;          % inlet conc of CO2 in gas
A       = 0.8855109225;  % area of the tray
k2      = 0.06;          % const of velocity rating
rhoMDEA = 1024.9;        % kg/m3
rhoH2O  = 992.21;        % kg/m3
MH2O    = 0.018;         % kg/mol
MMDEA   = 0.11916;       % kg/mol
xMDEA   = 0.1313868613;  % inlet conc of MDEA in liq
N       = 6;             % number of stages

pars = [T p wMDEA u S nLin nVin xCO2in yCO2in A k2 rhoMDEA rhoH2O MH2O MMDEA xMDEA N];

yguess = 0.5;
yguess = ones(N*7,1)*yguess;

ysol = mmfsolve(@(y) abs_react(y,pars), yguess); 

% 1col=nCO2; 2col=xCO2b; 3col=xCO2s; 4col=yCO2s; 5col=yCO2b; 6col=nV; 7col=nL

% converting calculating vector to matrix
NN=7; % number of unknown variables
ysoll = v2a(ysol,N,NN); disp(ysoll);

nCO2  = ysoll(:,1); yCO2b = ysoll(:,2); yCO2s = ysoll(:,3);
xCO2s = ysoll(:,4); xCO2b = ysoll(:,5); nV    = ysoll(:,6); nL = ysoll(:,7);

aa = zeros(N,1);
j=1;
for i=1:N
    aa(i) = j;
    j=j+1;
end


figure(1);
subplot(2,2,1);
plot(nCO2,aa,'g'); title('Molar Flow of CO2');
xlabel('nCO2 [mol/s]'); ylabel('Stages [-]');

subplot(2,2,3);
plot(nV,aa,'k'); title('Gas Flow [mol/s]');
xlabel('Gas Flow [mol/s]'); ylabel('Stages [-]');

subplot(2,2,4);
plot(nL,aa,'b'); title('Liquid Flow [mol/s]');
xlabel('Liquid Flow [mol/s]'); ylabel('Stages [-]');

figure(2);
subplot(2,2,1);
plot(yCO2b,aa,'m'); title('Mol. Frac. of CO2 in Bulk of Gas Phase');
xlabel('yCO2b [-]'); ylabel('Stages [-]');

subplot(2,2,2);
plot(yCO2s,aa,'c'); title('Mol. Frac. of CO2 on Surface in Gas Phase');
xlabel('yCO2s [-]'); ylabel('Stages [-]');

subplot(2,2,3);
plot(xCO2s,aa,'r'); title('Mol. Frac. of CO2 on Surface in Liquid Phase');
xlabel('xCO2s [-]'); ylabel('Stages [-]');

subplot(2,2,4);
plot(xCO2b,aa,'g'); title('Mol. Frac. of CO2 in Bulk of Liquid Phase');
xlabel('xCO2b [-]'); ylabel('Stages [-]');

