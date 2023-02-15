
function [KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES]=processtableau(Tableau,pH,pe)

n=size(Tableau,2);
A=cell2mat(Tableau(:,1:n-3));
K=cell2mat(Tableau(:,n-2));
P=cell2mat(Tableau(:,n-1));
SPECIESNAMES=strvcat(Tableau(:,n));% size(SPECIESNAMES)

C1=0; C2=0;
for i=1:size(P,1)
    if P(i)==0 
        C1=C1+1; 
        ASOLUTION(C1,:)=A(i,:);
        KSOLUTION(C1,:)=K(i,:);
        SOLUTIONNAMES(C1,:)=SPECIESNAMES(i,:);
    end
    if P(i)==1 
        C2=C2+1; 
        ASOLID(C2,:)=A(i,:);
        KSOLID(C2,:)=K(i,:);
        SOLIDNAMES(C2,:)=SPECIESNAMES(i,:);
    end
end

tstpH=isnan(pH); tstpe=isnan(pe);

oKSOLUTION=KSOLUTION;
oKSOLID=KSOLID;
oASOLUTION=ASOLUTION;
oASOLID=ASOLID;

save originaltableau.mat oKSOLUTION oKSOLID oASOLUTION oASOLID

if tstpH==0
    % adjust for fixed pH
[KSOLUTION,KSOLID,ASOLUTION,ASOLID]=get_equilib_fixed_pH(KSOLUTION,KSOLID,ASOLUTION,ASOLID,pH);
end

if tstpe==0
% adjust for fixed pe
[KSOLUTION,KSOLID,ASOLUTION,ASOLID]=get_equilib_fixed_pe(KSOLUTION,KSOLID,ASOLUTION,ASOLID,pe);
end


end