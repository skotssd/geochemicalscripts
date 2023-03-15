function [Cu,CuOH,CuOH2s,CuOs,tenorite,MASSERR]=CuOHtableau(pH,pe,T,flag1,flag2,flag3,flag4,flag5,database,IS)

% determine K values based on NIST and measured values versus ionic strength

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

Tableau=[...
%H e  Cu2+  logK         phase  species 
1  0  0     0            0      {'H'}
0  1  0     0            0      {'e'}
0  0  1     0            0      {'Cu'}
%end of identity matrix part
-1 0  0     logKw        0      {'OH'}
-1 0  1     logKh1       0      {'CuOH'}
%solids
-2 0  1     logKfCuOH2s  1      {'CuOH2s'}
-2 0  1     logKfCuOs    1      {'CuOs'}
-2 0  1     logKtenorite  1     {'tenorite'}
];

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5,database);

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end
%CuOs=0;
Cusolids=CuOH2s+CuOs+tenorite;  MASSERR=max(MASSERR);

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