function [res] = abs_react(yvec,pars)
% ABS_REACT
% CH4+CO2 (g), MDEA+H2O (l), transfer of CO2 with reaction in liquid bulk

T       = pars(1);  % K
p       = pars(2);  % Pa
wMDEA   = pars(3);  % weight fract of MDEA in liq
u       = pars(4);  % superficial velo. m/s
S       = pars(5);  % submergence of the liquid m
nLin    = pars(6);  % inlet liquid flow mol/s
nVin    = pars(7);  % inlet gas flow mol/s
xCO2in  = pars(8);  % inlet conc of CO2 in liq
yCO2in  = pars(9);  % inlet conc of CO2 in gas
A       = pars(10); % area of the tray
k2      = pars(11); % const of velocity rating
rhoMDEA = pars(12); % kg/m3
rhoH2O  = pars(13); % kg/m3
MH2O    = pars(14); % kg/mol
MMDEA   = pars(15); % kg/mol
xMDEA   = pars(16); % inlet conc of MDEA in liq
N       = pars(17); % number of MEZIPATER stages

j=1;
for i=1:N
    nCO2(i)  = yvec(j+0);
    yCO2b(i) = yvec(j+1);
    yCO2s(i) = yvec(j+2);
    xCO2s(i) = yvec(j+3);
    xCO2b(i) = yvec(j+4);
    nV(i)    = yvec(j+5);
    nL(i)    = yvec(j+6);
    j = j+7;
end


% mass transfer in the gas phase
nyCO2 = 26.9*10^(-6);                                                % molar volume m3/mol
nyCH4 = (15.9+2.31*4)*10^(-6);
MCO2  = (12+2*16)*10^(-3);                                           % molar weight kg/mol
MCH4  = (12+4*1)*10^(-3);
Dg    = 3.2*10^(-8)*T^(1.75) / (p*(nyCO2^(1/3)+nyCH4^(1/3))^2) * ...
        (1/MCO2 + 1/MCH4)^(0.5);                                     % m2/s
kg    = 0.000467 * (u*100)^(0.25) * (S*100)^(-0.5) * ...
        (Dg*10000)^(0.5);                                            % mol/cm2/s/atm

    
% mass transfer in the liquid phase
Dl = 6*10^(-13)*T^2 - 3*10^(-10)*T + 4*10^(-8);               % m2/s
kl = 13 * (u*100)^(0.25) * (S*100)^(-0.5) * (Dl*10000)^(0.5); % cm/s


% equilibrium
AA = -8.55445; BB = 4.01195; CC = 9.52345;       % param for Henry's law
pH2O  = 7386.34;                                 % vapour pressure Pa
TC    = 647;                                     % critic temp of solvent = H2O
TR    = T/TC;
tau   = 1-TR;
kHCO2 = pH2O*exp(AA/TR + BB*tau^(0.355)/TR + ...
        CC*TR^(-0.41)*exp(tau));                 % Henry Pa
KCO2  = kHCO2/p;                                 % distributing coef Pa/Pa


% effective interficial area
a = 0.535 * (u*100)^(0.25) * (S*100)^(0.83); % cm2/cm2
a = a*A;                                     % m2


% total conc of liquid
rhol = wMDEA*rhoMDEA + (1-wMDEA)*rhoH2O;      % density of liq kg/m3
Ml   = xMDEA*MMDEA + (1-xMDEA)*MH2O;          % molar weight of liq kg/mol
cTOT = rhol * ((1-wMDEA)/MH2O + wMDEA/MMDEA); % mol/m3

% inlet conditions
nL0    = nLin;
nV0    = nVin;
xCO2b0 = xCO2in;
yCO2b0 = yCO2in;


j=1;
% first stage
res(j+0) = nL0*xCO2b0 + nCO2(1) - k2*cTOT*xCO2b(1)*S*A - (nL(1)*xCO2b(1)); % CO2 bal in liq
res(j+1) = nV(2)*yCO2b(2) - (nV(1)*yCO2b(1) + nCO2(1));                    % CO2 bal in gas
res(j+2) = nCO2(1) - (kl*rhol*10^(-6)/Ml * (xCO2s(1) - ...
           xCO2b(1))*a*10^(4));                                            % flow CO2 balance
res(j+3) = kl*rhol*10^(-6)/Ml * (xCO2s(1) - xCO2b(1)) - ...
           kg*p*10^(-5)*(yCO2b(1) - yCO2s(1));                             % flow of CO2 through phases
res(j+4) = yCO2s(1) - (KCO2*xCO2s(1));                                     % eq on phase boundary
res(j+5) = nL(1) - (nL0 + nCO2(1));                                        % liq bal
res(j+6) = nV(2) - (nCO2(1) + nV(1));                                      % gas bal

j=8; % number of place in res vector
for i=2:N-1
    res(j+0) = nL(i-1)*xCO2b(i-1) + nCO2(i) - k2*cTOT*xCO2b(i)*S*A - ...
               (nL(i)*xCO2b(i));                                           % CO2 bal in liq
    res(j+1) = nV(i+1)*yCO2b(i+1) - (nCO2(i) + nV(i)*yCO2b(i));            % CO2 bal in gas
    res(j+2) = nCO2(i) - (kl*rhol*10^(-6)/Ml * (xCO2s(i) - ...
               xCO2b(i))*a*10^(4));                                        % flow CO2 bal
    res(j+3) = kl*rhol*10^(-6)/Ml * (xCO2s(i) - xCO2b(i)) - ...
               (kg*p*10^(-5) * (yCO2b(i) - yCO2s(i)));                     % flow of CO2 through phases
    res(j+4) = yCO2s(i) - (KCO2*xCO2s(i));                                 % eq on phase boundary
    res(j+5) = nL(i) - (nL(i-1) + nCO2(i));                                % liq bal
    res(j+6) = nV(i+1) - (nCO2(i) + nV(i));                                % gas bal
    j = j+7;
end

j = N*7-6; % number of place in res vector
% bottom stage (the last one)
res(j+0) = nL(N-1)*xCO2b(N-1) + nCO2(N) - k2*cTOT*xCO2b(N)*S*A - ...
           nL(N)*xCO2b(N);
res(j+1) = nV0*yCO2b0 - (nCO2(N) + nV(N)*yCO2b(N));
res(j+2) = nCO2(N) - (kl*rhol*10^(-6)/Ml * (xCO2s(N) - ...
           xCO2b(N))*a*10^(4));
res(j+3) = kl*rhol*10^(-6)/Ml * (xCO2s(N) - xCO2b(N)) - ...
           kg*p*10^(-5) * (yCO2b(N) - yCO2s(N));
res(j+4) = yCO2s(N) - KCO2*xCO2s(N);
res(j+5) = nL(N) - (nL(N-1) + nCO2(N));
res(j+6) = nV0 - (nCO2(N) + nV(N));

res = res';

end
