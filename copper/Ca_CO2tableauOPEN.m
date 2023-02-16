function [H2CO3,HCO3,CO3,CaCO3s,MASSERR]=Ca_CO2tableauOPEN(pH,pe,PCO2,T,flag1,flag2,flag3,flag4,flag5,database)

% input tableau.  change this part % ----------------------------------------------

logKH=-1.5; logKa1=-6.3; logKa2=-10.3; logKsp=-8.5;

Tableau=[...
%H      e        CO2g      Ca    logK                         phase    species1 
1       0        0         0      0                           0    {'H'}
0       1        0         0      0                           0    {'e '}
0       0        1         0      0                           -1   {'CO2g'}
0       0        0         1      0                           0    {'Ca'}
-1      0        0         0      -14                         0    {'OH'}
0       0        1         0      logKH                       0    {'H2CO3'}
-1      0        1         0      logKH+logKa1                0    {'HCO3'}
-2      0        1         0      logKH+logKa2+logKa1         0    {'CO3'}
-2      0        1         1      logKH+logKa2+logKa1-logKsp  1    {'CaCO3s'}
];

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableauOPEN(Tableau,pH,pe,PCO2);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5,database);

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

end