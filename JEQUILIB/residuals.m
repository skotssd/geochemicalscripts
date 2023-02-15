
function R=residuals(X,Asolution,Ksolution,Asolid,Ksolid,T)

Nx=size(Asolution,2); Ncp=size(Asolid,1); Nc=size(Asolution,1); 
Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);

% mass balance with only positive Xsolid values
Xsolidzero=Xsolid;
Xsolidzero(Xsolidzero < 0) = 0;
logC=Ksolution+Asolution*log10(Xsolution); C=10.^(logC); % calc species
Rmass=Asolution'*C+Asolid'*Xsolidzero-T;

% two versions of RSI
Q=Asolid*log10(Xsolution); SI=(Q+Ksolid);
for i=1:size(Xsolid,1)
    if Xsolid(i)>0; RSI(i)=(SI(i)); end % this should be zero if solids present
    if Xsolid(i)<=0 
        RSI(i)=(SI(i))-Xsolid(i); 
    end
end

R=[Rmass; RSI']; 

end