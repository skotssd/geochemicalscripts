function [Cu,CuOH,CuOH2s,CuOs,CuCO3s,tenorite,malachite,MASSERR]=CuOHCO2opentableau(pH,pe,PCO2,T,flag1,flag2,flag3,flag4,flag5,database)

% input tableau.  change this part % ----------------------------------------------

logKH=-1.5; logKa1=-6.3; logKa2=-10.3; logKsp=-8.5;

% determine K values based on NIST and measured values versus ionic strength
IS=0.01;  % eventually expand so this is a variable, and all values are interpolated.
% rxn Kw
logKw=H2O(IS);

% rxn Cu+OH=CuOH
logKfOH1=CuOH(IS);
logKh1=logKfOH1+logKw;

% solid phases have limited IS data in NIST

logKspCuOH2=-18.7; %zero ionic strength value
logKspCuO=-19.5; %zero ionic strength value
logKfCuOH2s=-1*logKspCuOH2+2*logKw;
logKfCuOs=-1*logKspCuO+2*logKw;
logKtenorite=20.18+2*logKw; %CuO
logKmalachite=33.18+2*logKw; %malachite Cu2CO3(OH)2
logKsp=-11.5; %CuCO3s

% Cu CO3 complexation

logKf1=6.47; logKf2=10.2;
logKfH1=1.03;

Tableau=[...
%H      e        CO2g      Cu    logK                                 phase    species1 
1       0        0         0      0                                   0    {'H'}
0       1        0         0      0                                   0    {'e '}
0       0        1         0      0                                  -1   {'CO2g'}
0       0        0         1      0                                   0    {'Cu'}
-1      0        0         0      logKw                               0    {'OH'}
0       0        1         0      logKH                               0    {'H2CO3'}
-1      0        1         0      logKH+logKa1                        0    {'HCO3'}
-2      0        1         0      logKH+logKa2+logKa1                 0    {'CO3'}
-1      0        0         1      logKh1                              0    {'CuOH'}
-2      0        1         1      logKH+logKa2+logKa1+logKf1          0    {'CuCO3'}
-4      0        2         1      2*logKH+2*logKa2+2*logKa1+logKf2    0    {'CuCO32'}
-1      0        1         1      logKH+logKa1+logKfH1                0    {'CuHCO3'}
%solids
-2      0        0         1      logKfCuOH2s                         1    {'CuOH2s'}
-2      0        0         1      logKfCuOs                           1    {'CuOs'}
-2      0        1         1      logKH+logKa2+logKa1-logKsp          1    {'CuCO3s'}
-2      0        0         1      logKtenorite                        1    {'tenorite'}
-2      0        1         2      logKH+logKa2+logKa1+logKmalachite   1    {'malachite'}
];

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableauOPEN(Tableau,pH,pe,PCO2);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5,database);

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

Cusolids=CuOH2s+CuOs+CuCO3s+tenorite+malachite;  MASSERR=max(MASSERR);

end

% ---------------- SUBFUNCTIONS --------------------------------------------------------

function logKh1=CuOH(ISvalue)

OH1=[...
0 6.5
0.1 6.1
0.5 6.1
0.7 6.2
1.0 6.3
2.0 6.6
3.0 (6.3+6.8)/2
3.00001 6.3
3.00002 6.8
];

IS=OH1(:,1); logK=OH1(:,2);
%output value
logKh1=interp1(IS,logK,ISvalue,'pchip');

end

function logKw=H2O(ISvalue)

IS=[0 0.1 0.5 0.5001 0.5002 0.7 1.0 1.00001 1.00002];
logK=-1*[13.997 13.78 13.73 13.69 13.75 13.75 13.77 13.71 13.94];

%output value
logKw=interp1(IS,logK,ISvalue,'pchip');

end