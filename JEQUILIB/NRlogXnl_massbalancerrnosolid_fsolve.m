
function [X,R,J,RSI,C] = NRlogXnl_massbalancerrnosolid_fsolve(Xguess,Asolution,Ksolution,Asolid,Ksolid,T,TYPX,flag2)

global Asolution Ksolution Asolid Ksolid T

options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',10000,...
    'TypicalX',TYPX,'Algorithm','trust-region-dogleg','StepTolerance',1e-25,'FunctionTolerance',1e-25,'MaxIterations',10000);
logXbest = fsolve(@(logX)speciation_logXnosolids(logX),log10(Xguess),options);

%function outputs
logX=logXbest; X=10.^logX;
[F,J]=speciation_logXnosolids(logXbest);
logC=(Ksolution)+Asolution*logX; C=10.^(logC); % calc species
R=Asolution'*C-T; 
Q=Asolid*log10(X); SI=(Q+Ksolid);
RSI=ones(size(SI))-SI;
F=[R]; 

end

%------------------- SUBFUNCTIONS ---------------------------------------

function [R,J]=speciation_logXnosolids(logX)

global Asolution Ksolution Asolid Ksolid T

[Nc,Nx]=size(Asolution); %Xsolution=X(1:Nx);
logC=(Ksolution)+Asolution*logX; C=10.^(logC); % calc species
R=Asolution'*C-T; 

% Evaluate the Jacobian of logX
% first the normal way like that CAryou paper
z=zeros(Nx,Nx);  X=10.^logX;
      
for j=1:Nx
	for k=1:Nx
		for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/X(k); end
   	end
end
% then multiply by X beause that is how you do the derivative wrt ln(X)
% and divide by 2.303 to convert to log base 10
for j=1:Nx
	for k=1:Nx
		z(j,k)=z(j,k)*(X(k)*2.303); 
   	end
end

J=z;

end
% 

% 
% [Nc,Nx]=size(Asolution); %Xsolution=X(1:Nx);
% criteria=1e-16; 
% 
% for II=1:100
% 
%   logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
%   R=Asolution'*C-T ;
% 
%    tester=isinf(X);
%    if max(tester)==1
%         for k=1:length(X)
%             tst=isinf(X(k));
%             if tst==1; X(k)=0.5*T(k); end
%         end
%         disp('no solid inf')
%         %X
%    end
% 
%    testthree=max(X./T); % for when mass balance is way crazy!
%    threshold=100;
%    if testthree>=threshold
%        %X
%        for k=1:Nx
%             if X(k)/T(k)>=threshold
%                 X(k)=0.5*T(k); 
%             end
%        end
%         disp('mass balance error logX no solids')
%    end
% 
%    testfour=max(abs(R)); % for when mass balance calculated R is way crazy!
%    threshold=1e15;
%    if testfour>=threshold
%        X=1e-35*X; disp('mass balance R error logX no solids')
%    end
% 
% 
% logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
% R=Asolution'*C-T ;
% 
% if flag2==1
% 
% % Evaluate the Jacobian of logX
% % first the normal way lilke that CAryou paper
%    z=zeros(Nx,Nx); 
%       
% for j=1:Nx
% 	for k=1:Nx
% 		for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/X(k); end
%    	end
% end
% % then multiply by X beause that is how you do the derivative wrt ln(X)
% % and divide by 2.303 to convert to log base 10
% 
% for j=1:Nx
% 	for k=1:Nx
% 		z(j,k)=z(j,k)*(X(k)*2.303); 
%    	end
% end
% 
% end
% 
% if flag2==2
% logX=log10(X);
% F=@(logX) residualslogXnosolids(logX,Asolution,Ksolution,Asolid,Ksolid,T);
% [z,err] = jacobianest(F,logX);
% %tst=max(isnan(z))
% end
% 
% logX=log10(X);
% % F=@(logX) residualslogXnosolids(logX,Asolution,Ksolution,Asolid,Ksolid,T);
% % [ztst,err] = jacobianest(F,logX);
% % z
% % ztst
% % pause
% 
% %z
% %log10(X)
% 
%     
%     %J=jac;
%     J=z ;
%     Jstar = J./TYPX';
%     %pause % scale Jacobian
%     %dx_star = -Jstar\(R); % faster than the alternatives it seems.  if warnings set to off.
%     %dx_star =pinv(z)*(-1*R); % no scaling
%     %dx_star =pinv(Jstar)*(-1*R);
%     %deltaX=J\(-1*R);
%     
%     %deltaX = pinv(J)*(-1*R);
%     %Jstar
%     dx_star = pinv(Jstar)*(-1*R);
%     %dx_star=Jstar\(-1*R);
%     deltaX = dx_star.*TYPX; %reverse the scaling
%     
%     %deltaX=dx;
%     % rescale x
% 
% %deltaX=z\(-1*R);
% %deltaX = pinv(z)*(-1*R);
% %one_over_del=max([1, -1*deltaX'./(0.5*X')]);
% %del=1/one_over_del; X=X+del*deltaX;
% logX=log10(X)+deltaX;
% X=10.^logX;
% %for loop=1:length(X(1:Nx)); if X(loop)>=T(loop); X(loop)=T(loop)/2; end; end
% 
% logC=(Ksolution)+Asolution*log10(X); C=10.^(logC); % calc species
% R=Asolution'*C-T ;
% 
% tst=sum(abs(R));
% 
% if tst<=criteria; break; end
% 
% end
% 
% Q=Asolid*log10(X); SI=(Q+Ksolid);
% RSI=ones(size(SI))-SI;
% 
% F=[R]; 
% 
% end