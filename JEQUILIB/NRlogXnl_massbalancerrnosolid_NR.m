
function [X,R,J,RSI,C] = NRlogXnl_massbalancerrnosolid_NR(X,Asolution,Ksolution,Asolid,Ksolid,T,TYPX,flag2,flag3)

Nx=size(Asolution,2); Ncp=size(Asolid,1); Nc=size(Asolution,1);  
Xsolution=X(1:Nx); %Xsolid=X(Nx+1:Nx+Ncp);
criteria=1e-16;  

for II=1:100
   %X
   logC=(Ksolution)+Asolution*log10(X);
   C=10.^(logC); % calc species
   R=Asolution'*C-T; R=R./T;
   
   tester=isinf(X);
   if max(tester)==1
        for k=1:length(X)
            tst=isinf(X(k));
            if tst==1; X(k)=0.99*T(k); end
        end
        if flag3==1; disp('no solid inf'); end
        %X
   end

   % for when zero in Xsolution.  can't take log of it.
   testtwo=ismember(0, Xsolution);
   
   if max(testtwo)==1
       for k=1:Nx
            if Xsolution(k)==0; Xsolution(k)=0.99*T(k); end
       end
       if flag3==1; disp('zero no solid'); end
   end

   testthree=max(X./T); % for when mass balance is way crazy!
   threshold=100;
   if testthree>=threshold
       %X
       for k=1:Nx
            if X(k)/T(k)>=threshold
                X(k)=0.99*T(k); 
            end
       end
        if flag3==1; disp('mass balance error logX no solids'); end
   end

   testfour=max(abs(R)); % for when mass balance calculated R is way crazy!
   threshold=1e15;
   if testfour>=threshold
       %X
       for k=1:Nx
            if X(k)/T(k)>=threshold
                %X(index)=0.99*T(index); 
                X(k)=0.99*T(k); 
                %X(k)=1e-20*X(k);
            end
       end
       if flag3==1; disp('mass balance R error logX no solids'); end
       %X
   end


logC=(Ksolution)+Asolution*log10(X); C=10.^(logC);
% calc species
R=Asolution'*C-T ;

if flag2==1

% Evaluate the Jacobian of logX
% first the normal way lilke that CAryou paper
   z=zeros(Nx,Nx); 
      
for j=1:Nx
	for k=1:Nx
		for i=1:Nc; z(j,k)=z(j,k)+(Asolution(i,j)*Asolution(i,k)*C(i)/X(k)); end
   	end
end
% then multiply by X beause that is how you do the derivative wrt ln(X)
% and divide by 2.303 to convert to log base 10
% z
for j=1:Nx
	for k=1:Nx
		z(j,k)=z(j,k)*(X(k)*log(10));
   	end
end

%     logX=log10(X)
%     F=@(logX) residualslogXnosolids(logX,Asolution,Ksolution,Asolid,Ksolid,T);
%     [zcompare,err] = jacobianest(F,logX);
%     z
%     zcompare
%     zcompare./z
%     pause 

end

if flag2==2
logX=log10(X);
F=@(logX) residualslogXnosolids(logX,Asolution,Ksolution,Asolid,Ksolid,T);
[z,err] = jacobianest(F,logX);
%tst=max(isnan(z))
end

logX=log10(X);
% F=@(logX) residualslogXnosolids(logX,Asolution,Ksolution,Asolid,Ksolid,T);
% [ztst,err] = jacobianest(F,logX);
% z
% ztst
% pause

%z
%log10(X)

    
    %J=jac;
    J=z ;
    Jstar = J./TYPX';
    %pause % scale Jacobian
    %dx_star = -Jstar\(R); % faster than the alternatives it seems.  if warnings set to off.
    %dx_star =pinv(z)*(-1*R); % no scaling
    %dx_star =pinv(Jstar)*(-1*R);
    %deltaX=J\(-1*R);
    
    %deltaX = pinv(J)*(-1*R);
    %Jstar
  
    dx_star = pinv(Jstar)*(-1*R);
    %dx_star=Jstar\(-1*R);
    deltaX = dx_star.*TYPX; %reverse the scaling
    
    %deltaX=dx;
    % rescale x

%deltaX=z\(-1*R);
%deltaX = pinv(z)*(-1*R);
%one_over_del=max([1, -1*deltaX'./(0.5*X')]);
%del=1/one_over_del; X=X+del*deltaX;
logX=log10(X)+deltaX;
X=10.^logX;
%for loop=1:length(X(1:Nx)); if X(loop)>=T(loop); X(loop)=T(loop)/2; end; end

logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
R=Asolution'*C-T ;

tst=sum(abs(R));

if tst<=criteria; break; end

end

Q=Asolid*log10(X); SI=(Q+Ksolid);
RSI=ones(size(SI))-SI;

F=[R]; 

end