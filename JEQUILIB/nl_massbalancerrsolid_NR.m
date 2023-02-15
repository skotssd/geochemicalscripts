
function [X,F,J,RSI,C] = nl_massbalancerrsolid_NR(X,Asolution,Ksolution,Asolid,Ksolid,T,TYPX,flag2)

Nx=size(Asolution,2); Ncp=size(Asolid,1); Nc=size(Asolution,1);  
Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);
criteria=1e-16; 

for i=1:100

    % mass balance with only positive Xsolid values
    Xsolidzero=Xsolid;
    Xsolidzero(Xsolidzero < 0) = 0; 
    logC=Ksolution+Asolution*log10(Xsolution); C=10.^(logC);
    % calc species
    Rmass=Asolution'*C+Asolid'*Xsolidzero-T;
    Rsolid=(Asolid'*Xsolidzero-T)./T;
    Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);

    tester=isinf(X); % for when X has an infinite value.  can't have that! so set to half of total mass
   
   if max(tester)==1
        for k=1:length(X)
            tst=isinf(X(k));
            if tst==1; X(k)=0.99*T(k); end
        end
        disp('inf')
   end
   
   % for when zero in Xsolution.  can't take log of it.
   testtwo=ismember(0, Xsolution);
   
   if max(testtwo)==1
       for k=1:Nx
            if Xsolution(k)==0; Xsolution(k)=0.99*T(k); end
       end
        disp('zero')
   end
    
   testthree=max(Xsolution./T); % for when mass balance is way crazy!
   threshold=1e2;
   if testthree>=threshold
       for k=1:Nx
            if Xsolution(k)/T(k)>=threshold
                Xsolution(k)=0.99*T(k); 
            end
       end
        disp('mass balance solidsversion')
   end
       
   testfour=max(abs(Rmass./T)); % for when mass balance calculated R is way crazy!
   threshold=1e15;
   if testfour>=threshold
       [value,index]=max(abs(Rmass./T));
         X(index)=0.99*T(index); 
   disp('mass balance R error logX solids')
   end


% mass balance with only positive Xsolid values
Xsolidzero=Xsolid;
Xsolidzero(Xsolidzero < 0) = 0;
logC=Ksolution+Asolution*log10(Xsolution); C=10.^(logC); % calc species
Rmass=Asolution'*C+Asolid'*Xsolidzero-T;
Q=Asolid*log10(Xsolution); SI=(Q+Ksolid);
RSI=SI;

% two versions of RSI
Q=Asolid*log10(Xsolution); SI=(Q+Ksolid); 
for i=1:size(Xsolid,1)
    if Xsolid(i)>0; RSI(i)=(SI(i)); end % this should be close to zero if solids present
    if Xsolid(i)<=0 
        RSI(i)=(SI(i))-Xsolid(i);
    end
end

R=[Rmass; RSI];

if flag2==1

% Evaluate the Jacobian 
   z=zeros(Nx+Ncp,Nx+Ncp);
	for j=1:Nx % mass balance error in terms of solution components
		for k=1:Nx 
				for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/Xsolution(k); end
       	end
    end
    for j=1:Nx % mass balance error terms by solid components
		for k=Nx+1:Nx+Ncp
            z(j,k)=Asolid(k-Nx,j); %derivative of mass balance error by each solid conc (when solid +ve)
            if Xsolid(k-Nx)<=0; z(j,k)=0; end %zero when not in mass balance.
       	end
    end
    for j=Nx+1:Nx+Ncp % RSI term by solution components
		for k=1:Nx
				%z(j,k)=-1*Asolid(j-Nx,k)*(SI(j-Nx)/Xsolution(k)); % for SI-1 formulation
                z(j,k)=Asolid(j-Nx,k)/(log(10)*Xsolution(k)); % for log(SI) formulation
                %z(j,k)=Asolid(j-Nx,k)/Xsolution(k); % for ln(SI)-Xsolid formulation
                %if Xsolid(k-Nx)<=0; z(j,k)=Asolid(j-Nx,k)/Xsolution(k); end; %same when not in mass balance.
      	end
    end
    for j=Nx+1:Nx+Ncp %RSI term by solid components
        for k=Nx+1:Nx+Ncp
            z(j,k)=0;
            if Xsolid(k-Nx)<=0; z(j,k)=0; if j==k; z(j,k)=-1; end; end % logSI-Xsolid when Xsolid is negative.  so derivative -1
        end
    end
  
end

if flag2==2
    F=@(X) residuals(X,Asolution,Ksolution,Asolid,Ksolid,T);
    [z,err] = jacobianest(F,X); z=real(z);
end

    J=z;
    
    %TYPX=ones(size(X));
    weight = ones(Nx+Ncp,1); 
    J0 = weight*(1./TYPX'); % Jacobian scaling matrix
    Jstar = J./J0;
    %pause % scale Jacobian
    %dx_star = -Jstar\(R); % might be faster but does not give good answers for solids
    %dx_star =pinv(z)*(-1*R); 
    dx_star =pinv(Jstar)*(-1*R); % seems to give the best answers
    %deltaX=z\(-1*R);
    %deltaX = pinv(Jstar)*(-1*R);
    %dx_star = pinv(Jstar)*(-1*R);
    dx = dx_star.*TYPX;
    deltaX=dx;
    % rescale x
 
%deltaX = pinv(z)*(-1*R);
Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);
delXsolution=deltaX(1:Nx); delXsolid=deltaX(Nx+1:Nx+Ncp);
one_over_del=max([1, -1*delXsolution'./(0.5*Xsolution')]);
del=1/one_over_del; %del=1;
Xsolution=Xsolution+del*delXsolution;
Xsolid=Xsolid+delXsolid;
X=[Xsolution; Xsolid];
         
tst=sum(abs(R));
if tst<=criteria; break; end
end

% Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);
% Q=Asolid*log10(Xsolution); SI=10.^(Q+Ksolid);
% RSI=ones(size(SI))-SI;
% 

% mass balance with only positive Xsolid values
Xsolidzero=Xsolid;
Xsolidzero(Xsolidzero < 0) = 0;
logC=Ksolution+Asolution*log10(Xsolution); C=10.^(logC); % calc species
Rmass=Asolution'*C+Asolid'*Xsolidzero-T;
Q=Asolid*log10(Xsolution); SI=(Q+Ksolid);
RSI=SI;

% two versions of RSI
Q=Asolid*log10(Xsolution); SI=(Q+Ksolid);
for i=1:size(Xsolid,1)
    if Xsolid(i)>0; RSI(i)=(SI(i)); end % this should be zero if solids present
    if Xsolid(i)<=0 
        %RSI(i)=0;
        RSI(i)=(SI(i))-Xsolid(i); 
        %RSI(i)=(SI(i))-Xsolid(i)-2; 
    end
end

F=[Rmass]; J=z;

end