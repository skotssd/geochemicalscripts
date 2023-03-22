
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,F,J,C] = nl_massbalancerrnosolid_NR(X,Asolution,Ksolution,T)

[Nc,Nx]=size(Asolution); %Xsolution=X(1:Nx);
criteria=1e-15;

for i=1:1000

logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
R=Asolution'*C-T;

% Evaluate the Jacobian 
   z=zeros(Nx,Nx); 
for j=1:Nx 
	for k=1:Nx 
		for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/X(k); end
   	end
end

J = z;

deltaX=z\(-1*R);
one_over_del=max([1, -1*deltaX'./(0.5*X')]);
del=1/one_over_del; X=X+del*deltaX;
    
tst=sum(abs(R));
if tst<=criteria; break; end

end

F=[R]; 

end