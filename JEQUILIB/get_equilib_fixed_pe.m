
% ----------- for fixed pe ----------------

function [Ksolution,Ksolid,Asolution,Asolid]=get_equilib_fixed_pe(KSOLUTION,KSOLID,ASOLUTION,ASOLID,pe)

    [N,M]=size(ASOLUTION);
    Ksolution=KSOLUTION-ASOLUTION(:,1)*pe;
    Asolution=[ASOLUTION(:,2:M)];
    [N,M]=size(ASOLID);
    Ksolid=KSOLID-ASOLID(:,1)*pe;
    Asolid=[ASOLID(:,2:M)];

end