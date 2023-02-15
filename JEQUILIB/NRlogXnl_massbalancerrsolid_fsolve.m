
function [X,F,J,RSI,C] = NRlogXnl_massbalancerrsolid_fsolve(Xguess,Asolution,Ksolution,Asolid,Ksolid,T,TYPX,flag2)

global Asolution Ksolution Asolid Ksolid T

Nx=size(Asolution,2); Ncp=size(Asolid,1); Nc=size(Asolution,1); 
Xtransform=[log10(Xguess(1:Nx)) Xguess(Nx+1:length(Xguess))]
[tst,Jtst]=speciation_logXsolids((Xtransform))
%pause
options = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',10000,...
    'TypicalX',TYPX,'Algorithm','trust-region-dogleg','StepTolerance',1e-25,'FunctionTolerance',1e-25,'MaxIterations',10000);
logXbest = fsolve(@(logX)speciation_logXsolids(logX),Xtransform,options);
[tst,Jtst]=speciation_logXsolids((logXbest))
%pause
logXbest = fsolve(@(logX)speciation_logXsolids(logX),logXbest,options);

%pause

Nx=size(Asolution,2); Ncp=size(Asolid,1); Nc=size(Asolution,1); 
Xsolution=10.^(logXbest(1:Nx)); Xsolid=logXbest(Nx+1:Nx+Ncp); X=logXbest;
 
% % mass balance with only positive Xsolid values
[F,J]=speciation_logXsolids((logXbest));
Xsolidzero=Xsolid;
Xsolidzero(Xsolidzero < 0) = 0; 
logC=Ksolution+Asolution*log10(Xsolution); C=10.^(logC);
% calc species
R=Asolution'*C+Asolid'*Xsolidzero-T; F=R;
Q=Asolid*log10(Xsolution); SI=(Q+Ksolid);
RSI=SI;

end

%------------------- SUBFUNCTIONS ---------------------------------------

function [R,J]=speciation_logXsolids(logX)

global Asolution Ksolution Asolid Ksolid T

Nx=size(Asolution,2); Ncp=size(Asolid,1); Nc=size(Asolution,1); 
Xsolution=10.^(logX(1:Nx)); Xsolid=logX(Nx+1:Nx+Ncp); X=logX;
 
% % mass balance with only positive Xsolid values
Xsolidzero=Xsolid;
Xsolidzero(Xsolidzero < 0) = 0; 
logC=Ksolution+Asolution*log10(Xsolution); C=10.^(logC);
% calc species
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

   
% % Evaluate the Jacobian on logX.  but first as normal
    z=zeros(Nx+Ncp,Nx+Ncp);
	for j=1:Nx % mass balance error in terms of solution components
		for k=1:Nx 
				for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/Xsolution(k); end
       	end
    end
    
    % for log ofX derivatives (just on the soluble part). the solution part
    % works.
    
    for j=1:Nx
	for k=1:Nx
        z(j,k)=z(j,k)*(X(k)*2.303); 
    end
    end
    %X
    %z
    for j=1:Nx % mass balance error terms by solid components
		for k=Nx+1:Nx+Ncp
            z(j,k)=Asolid(k-Nx,j); %derivative of mass balance error by each solid conc (when solid +ve)
            if Xsolid(k-Nx)<=0; z(j,k)=0; end %zero when not in mass balance.
       	end
    end
    
       
    for j=Nx+1:Nx+Ncp % RSI term by solution components (no change for log of X)
		for k=1:Nx
				%z(j,k)=-1*Asolid(j-Nx,k)*(SI(j-Nx)/Xsolution(k)); % for SI-1 formulation
                z(j,k)=((Asolid(j-Nx,k)/(log(10)*Xsolution(k))))*(X(k)*2.303); % for log(SI) formulation
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

J=z;

end


% Nx=size(Asolution,2); Ncp=size(Asolid,1); Nc=size(Asolution,1);  %TYPX(1:Nx)=(log10(TYPX(1:Nx)));
% 
% typX=ones(size(X)); typX(1:Nx)=((TYPX(1:Nx))); TYPX=typX;
% 
% Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);
% 
% %Rtst=residuals(X,Asolution,Ksolution,Asolid,Ksolid,T);
% %fun = @(X) residuals(X,Asolution,Ksolution,Asolid,Ksolid,T);
% %[jac,err] = jacobianest(fun,X);
% %jac
% %pause
% criteria=1e-16; 
% 
% %for i=1:100
% tst=1; %initialize variable for terminate loop
% COUNT=0;
% 
% while tst>=criteria
%     COUNT=COUNT+1;
%     
%    tester=isinf(X); % for when X has an infinite value.  can't have that! so set to half of total mass
%    
%    if max(tester)==1
%         for k=1:length(X)
%             tst=isinf(X(k));
%             if tst==1; X(k)=0.5*T(k); end
%         end
%         disp('inf')
%         X
%    end
%    
%    % for when X is greater than T.  can't have that! so set to half of total mass
%    %tester=Xsolution-T(1:length(Xsolution));
%       
%    %if max(tester)>0
%    %     for k=1:length(X)
%    %         tst=isinf(X(k));
%    %         if tst==1; X(k)=0.5*T(k); end
%    %     end
%    %     disp('mass balance error')
%    %     X
%    %end
%            
%    Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);
% 
%    testtwo=Xsolution(Xsolution==0); % for when zero in Xsolution.  can't take log of it.
%    
%    if testtwo==0
%        for k=1:Nx
%             if Xsolution(k)==0; Xsolution(k)=1e-10*T(k); end
%        end
%         disp('zero')
%         Xsolution
%    end
%     
%    testthree=max(Xsolution./T); % for when mass balance is way crazy!
%    threshold=10;
%    if testthree>=threshold
%        Xsolution
%        for k=1:Nx
%             if Xsolution(k)/T(k)>=threshold; 
%                 Xsolution(k)=0.5*T(k); 
%             end
%             %Xsolid=-6*ones(size(Xsolid));
%        end
%         disp('mass balance solidsversion')
%         Xsolution
%    end
%        
%    
%    
% % mass balance with only positive Xsolid values
% Xsolidzero=Xsolid;
% Xsolidzero(Xsolidzero < 0) = 0; 
% logC=Ksolution+Asolution*log10(Xsolution); C=10.^(logC);
% % calc species
% Rmass=Asolution'*C+Asolid'*Xsolidzero-T;
% Q=Asolid*log10(Xsolution); SI=(Q+Ksolid);
% RSI=SI;
% 
% % two versions of RSI
% Q=Asolid*log10(Xsolution); SI=(Q+Ksolid); 
% for i=1:size(Xsolid,1)
%     if Xsolid(i)>0; RSI(i)=(SI(i)); end % this should be close to zero if solids present
%     if Xsolid(i)<=0 
%         RSI(i)=(SI(i))-Xsolid(i);
%     end
% end
% 
% R=[Rmass; RSI];
% tst=sum(abs(R));
% 
% if COUNT>=100; tst=0; disp('logX solid ITER EXEED'); end % just  make sure don't  have infinite loop
% 
% %if tst<=criteria; disp('should end'); pause; end
% 
% 
% if flag2==1
% 
% % Evaluate the Jacobian on logX.  but first as normal
%    z=zeros(Nx+Ncp,Nx+Ncp);
% 	for j=1:Nx % mass balance error in terms of solution components
% 		for k=1:Nx 
% 				for i=1:Nc; z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/Xsolution(k); end
%        	end
%     end
%     
%         % for log ofX derivatives (just on the soluble part). the solution part
%     % works.
%     
%     for j=1:Nx
% 	for k=1:Nx
%         z(j,k)=z(j,k)*(X(k)*2.303); 
%     end
%     end
%     %X
%     %z
%     for j=1:Nx % mass balance error terms by solid components
% 		for k=Nx+1:Nx+Ncp
%             z(j,k)=Asolid(k-Nx,j); %derivative of mass balance error by each solid conc (when solid +ve)
%             if Xsolid(k-Nx)<=0; z(j,k)=0; end %zero when not in mass balance.
%        	end
%     end
%     
%        
%     for j=Nx+1:Nx+Ncp % RSI term by solution components (no change for log of X)
% 		for k=1:Nx
% 				%z(j,k)=-1*Asolid(j-Nx,k)*(SI(j-Nx)/Xsolution(k)); % for SI-1 formulation
%                 z(j,k)=((Asolid(j-Nx,k)/(log(10)*Xsolution(k))))*(X(k)*2.303); % for log(SI) formulation
%                 %z(j,k)=Asolid(j-Nx,k)/Xsolution(k); % for ln(SI)-Xsolid formulation
%                 %if Xsolid(k-Nx)<=0; z(j,k)=Asolid(j-Nx,k)/Xsolution(k); end; %same when not in mass balance.
%       	end
%     end
%     for j=Nx+1:Nx+Ncp %RSI term by solid components
%         for k=Nx+1:Nx+Ncp
%             z(j,k)=0;
%             if Xsolid(k-Nx)<=0; z(j,k)=0; if j==k; z(j,k)=-1; end; end % logSI-Xsolid when Xsolid is negative.  so derivative -1
%         end
%     end
%     
%     
% end
% 
%     logX=[log10(Xsolution)
%         Xsolid];
%     
% if flag2==2
%     F=@(logX) residualslogXwithsolids(logX,Asolution,Ksolution,Asolid,Ksolid,T);
% 
%     [z,err] = jacobianest(F,logX);
% end
% 
%     %F=@(logX) residualslogXwithsolids(logX,Asolution,Ksolution,Asolid,Ksolid,T);
% 
%    % [ztst,err] = jacobianest(F,logX);
% 
% %     jac
%     % z
%     % ztst
%     % pause
% %     logX
%     %pause
%     %check each term of jacobian
%     %one=[z(1,1) jac(1,1)]
%     %two=[z(1,2) jac(1,2)]
%     %three=[z(2,1) jac(2,1)]
%     %four=[z(2,2) jac(2,2)]
%     
%     %pause
%     
%     %J=jac;
%     J=z;
%     %TYPX=logX;
%     
%     Jstar = J./TYPX';
%     
%     %pause % scale Jacobian
%     %dx_star = -Jstar\(R); % faster than the alternatives it seems.  if warnings set to off.
%     %dx_star =pinv(z)*(-1*R); % no scaling
%     %dx_star =pinv(Jstar)*(-1*R);
%     %deltaX=J\(-1*R);
%     
%     %without scaling ------------
%   
% %    deltaX = pinv(J)*(-1*R);
%     
%     % with scaling ---------------
%      dx_star = pinv(Jstar)*(-1*R);
%      %dx_star = -Jstar\(R);
%      %TEST=isinf(dx_star);
%      %if max(TEST)==1; dx_star=dx_star1; end
%      dx = dx_star.*TYPX; %reverse the scaling
%      deltaX=dx;
%     % rescale x
% 
% %deltaX=z\(-1*R);
% %deltaX = pinv(z)*(-1*R);
% %one_over_del=max([1, -1*deltaX'./(0.5*X')]);
% %del=1/one_over_del; X=X+del*deltaX;
% %Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);
% 
% logX=logX+deltaX; 
% X=[10.^logX(1:Nx) 
%     logX(Nx+1:Nx+Ncp)];
% %for loop=1:length(X(1:Nx)); if X(loop)>=T(loop); X(loop)=T(loop)/2; end; end
% %X
% %pause
% 
% 
% end
% 
% % Xsolution=X(1:Nx); Xsolid=X(Nx+1:Nx+Ncp);
% % Q=Asolid*log10(Xsolution); SI=10.^(Q+Ksolid);
% % RSI=ones(size(SI))-SI;
% % 
% 
% % mass balance with only positive Xsolid values
% Xsolidzero=Xsolid;
% Xsolidzero(Xsolidzero < 0) = 0;
% logC=Ksolution+Asolution*log10(Xsolution); C=10.^(logC); % calc species
% Rmass=Asolution'*C+Asolid'*Xsolidzero-T;
% Q=Asolid*log10(Xsolution); SI=(Q+Ksolid);
% RSI=SI;
% 
% % two versions of RSI
% Q=Asolid*log10(Xsolution); SI=(Q+Ksolid);
% for i=1:size(Xsolid,1)
%     if Xsolid(i)>0; RSI(i)=(SI(i)); end % this should be zero if solids present
%     if Xsolid(i)<=0 
%         %RSI(i)=0;
%         RSI(i)=(SI(i))-Xsolid(i); 
%         %RSI(i)=(SI(i))-Xsolid(i)-2; 
%     end
% end
% 
% F=[Rmass]; J=z;
% 
% end