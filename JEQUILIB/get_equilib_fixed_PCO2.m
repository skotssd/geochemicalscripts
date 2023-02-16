
% ----------- for fixed PCO2 ----------------

function [Ksolution,Ksolid,Asolution,Asolid]=get_equilib_fixed_PCO2(KSOLUTION,KSOLID,ASOLUTION,ASOLID,PCO2)

    [N,M]=size(ASOLUTION);
    Ksolution=KSOLUTION+ASOLUTION(:,1)*log10(PCO2);
    Asolution=[ASOLUTION(:,2:M)];
    [N,M]=size(ASOLID);
    Ksolid=KSOLID+ASOLID(:,1)*log10(PCO2);
    Asolid=[ASOLID(:,2:M)];

end