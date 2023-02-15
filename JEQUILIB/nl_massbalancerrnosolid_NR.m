

function [X,F,J,RSI,C] = nl_massbalancerrnosolid_NR(X,Asolution,Ksolution,Asolid,Ksolid,T,TYPX)

[Nc,Nx]=size(Asolution); %Xsolution=X(1:Nx);

criteria=1e-16;

tst=1; %initialize variable for terminate loop
COUNT=0;

while tst>=criteria
    COUNT=COUNT+1;
logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
R=Asolution'*C-T; 
    tester=isinf(X);
   if max(tester)==1
        for k=1:length(X)
            tst=isinf(X(k));
            if tst==1; X(k)=0.5*T(k); end
        end
        disp('no solid inf')
        %X
   end

   testthree=max(X./T); % for when mass balance is way crazy!
   threshold=10;
   if testthree>=threshold
       %X
       for k=1:Nx
            if X(k)/T(k)>=threshold
                X(k)=0.5*T(k); 
            end
       end
        disp('mass balance X no solids')
        
   end
       
   testfour=max(abs(R)); % for when mass balance calculated R is way crazy!
   threshold=1e2;
   if testfour>=threshold
       %X=1e-35*X; disp('mass balance R error logX no solids')
       for k=1:Nx
            if X(k)/T(k)>=threshold
                X(k)=0.01*T(k); 
            end
       end
       disp('mass balance R error X no solids')
   end

   


logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
R=Asolution'*C-T; 

% Evaluate the Jacobian 
   z=zeros(Nx,Nx); 
for j=1:Nx; 
	for k=1:Nx; 
		for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/X(k); end
   	end
end

    %[jac,err] = jacobianest(fun,X);
    %J=jac;
    J=z;
    Jstar = J./TYPX';
    %pause % scale Jacobian
    dx_star = -Jstar\(R); % faster than the alternatives it seems.  if warnings set to off.
    %dx_star =pinv(z)*(-1*R); % no scaling
    %dx_star =pinv(Jstar)*(-1*R);
    %deltaX=z\(-1*R);
    %deltaX = pinv(Jstar)*(-1*R);
    %dx_star = pinv(Jstar)*(-1*R);
    dx = dx_star.*TYPX; %reverse the scaling
    deltaX=dx;
    % rescale x

%deltaX=z\(-1*R);
%deltaX = pinv(z)*(-1*R);
one_over_del=max([1, -1*deltaX'./(0.5*X')]);
del=1/one_over_del; X=X+del*deltaX;
    
tst=sum(abs(R));

if COUNT>=100; tst=0; disp('X no ITER EXEED no solid logX'); end % just  make sure don't  have infinite loop

end

Q=Asolid*log10(X); SI=(Q+Ksolid);
RSI=ones(size(SI))-SI;

F=[R]; 

end