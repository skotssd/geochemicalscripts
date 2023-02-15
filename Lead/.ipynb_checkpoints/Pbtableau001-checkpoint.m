function [Pb,PbOH,PbCl,Pbsolids,MASSERR]=Pbtableau001(pH,pe,T,flag1,flag2,flag3,flag4,flag5,database)

Tableau=[...
%H,e,Pb,Cl,CO3,Ac,Y,logK,phase,species 
1,0,0,0,0,0,0,0,0,{'H'}
0,1,0,0,0,0,0,0,0,{'e'}
0,0,1,0,0,0,0,0,0,{'Pb'}
0,0,0,1,0,0,0,0,0,{'Cl'}
0,0,0,0,1,0,0,0,0,{'CO3'}
0,0,0,0,0,1,0,0,0,{'Ac'}
0,0,0,0,0,0,1,0,0,{'Y'}
%end of identity matrix
-1,0,0,0,0,0,0,-13.9944,0,{'OH'}
-1,0,1,0,0,0,0,-7.6005,0,{'PbOH'}
-2,0,1,0,0,0,0,-17.0898,0,{'PbOH2'}
-3,0,1,0,0,0,0,-28.0832,0,{'PbOH3'}
-2,0,1,0,0,0,0,-12.6888,1,{'PbOH2s'}
0,0,1,1,0,0,0,1.5546,0,{'PbCl'}
0,0,1,2,0,0,0,1.8953,0,{'PbCl2'}
0,0,1,3,0,0,0,1.7947,0,{'PbCl3'}
0,0,1,2,0,0,0,4.7801,1,{'PbCl2s'}
1,0,0,0,1,0,0,10.3239,0,{'HCO3'}
2,0,0,0,1,0,0,16.6743,0,{'H2CO3'}
0,0,1,0,1,0,0,5.4,0,{'PbCO3'}
0,0,1,0,1,0,0,13.1847,1,{'cerrusite'}
0,0,1,0,2,0,0,8.86,0,{'PbCO32'}
1,0,1,0,1,0,0,12.2339,0,{'PbHCO3'}
-1,0,1,0,1,0,0,-3.0944,0,{'PbOHCO3'}
-2,0,3,0,2,0,0,15.8112,1,{'hydrocerrusite'}
1,0,0,0,0,1,0,4.7542,0,{'Hac'}
0,0,1,0,0,1,0,2.5747,0,{'PbAc'}
0,0,1,0,0,2,0,4.012,0,{'PbAc2'}
0,0,1,0,0,3,0,3.4652,0,{'PbAc3'}
1,0,0,0,0,0,1,9.822,0,{'HY'}
2,0,0,0,0,0,1,16.0933,0,{'H2Y'}
3,0,0,0,0,0,1,18.7831,0,{'H3Y'}
4,0,0,0,0,0,1,20.7851,0,{'H4Y'}
5,0,0,0,0,0,1,22.3046,0,{'H5Y'}
6,0,0,0,0,0,1,22.1048,0,{'H6Y'}
0,0,1,0,0,0,1,18,0,{'PbY'}
1,0,1,0,0,0,1,12.622,0,{'PbHY'}
2,0,1,0,0,0,1,17.7933,0,{'PbH2Y'}
3,0,1,0,0,0,1,19.9831,0,{'PbH3Y'}
];

% end of tableau.  ------------------ % ----------------------------------------------

[KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe);

[SPECIESCONCS,SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5,database);

for k=1:size(SPECIESCONCS,1)
      txt=[SPECIATIONNAMES(k,:),'=SPECIESCONCS(k);'];
      eval(txt)
end

Pbsolids=PbOH2s+PbCl2s+cerrusite+hydrocerrusite;  MASSERR=max(MASSERR);

end