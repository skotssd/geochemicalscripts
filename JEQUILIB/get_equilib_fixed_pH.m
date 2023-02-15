
function [Ksolution,Ksolid,Asolution,Asolid]=get_equilib_fixed_pH(KSOLUTION,KSOLID,ASOLUTION,ASOLID,pH)

    [N,M]=size(ASOLUTION);
    Ksolution=KSOLUTION-ASOLUTION(:,1)*pH;
    Asolution=[ASOLUTION(:,2:M)];
    [N,M]=size(ASOLID);
    Ksolid=KSOLID-ASOLID(:,1)*pH;
    Asolid=[ASOLID(:,2:M)];

end