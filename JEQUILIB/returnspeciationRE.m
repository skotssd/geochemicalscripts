
function [C, SPECIATIONNAMES,MASSERR,X]=returnspeciationRE(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5,database,acid)

% NR on X with either analytical or numerical gradients --------------

if flag1==1
     for i=1:size(SOLUTIONNAMES,1)
         tst=SOLUTIONNAMES(i,:);
         for j=1:length(tst)
             if tst(j)=='+'; tst(j)='p'; end
             if tst(j)=='-'; tst(j)='m'; end
             if tst(j)=='('; tst(j)='L'; end
             if tst(j)==')'; tst(j)='R'; end
             if tst(j)==':'; tst(j)='C'; end
         end
         SOLUTIONNAMES(i,:)=tst;
     end
     
     for i=1:size(SOLIDNAMES,1)
         tst=SOLIDNAMES(i,:);
         for j=1:length(tst)
             if tst(j)=='('; tst(j)='L'; end
             if tst(j)==')'; tst(j)='R'; end
             if tst(j)==':'; tst(j)='C'; end
         end
         SOLIDNAMES(i,:)=tst;
     end
     
     [C, SPECIATIONNAMES,MASSERR,X]=NRX(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2);
end

if flag1==2
    
    for i=1:size(SOLUTIONNAMES,1)
         tst=SOLUTIONNAMES(i,:);
         for j=1:length(tst)
             if tst(j)=='+'; tst(j)='p'; end
             if tst(j)=='-'; tst(j)='m'; end
             if tst(j)=='('; tst(j)='L'; end
             if tst(j)==')'; tst(j)='R'; end
         end
         SOLUTIONNAMES(i,:)=tst;
     end
     
     for i=1:size(SOLIDNAMES,1)
         tst=SOLIDNAMES(i,:);
         for j=1:length(tst)
             if tst(j)=='('; tst(j)='L'; end
             if tst(j)==')'; tst(j)='R'; end
             if tst(j)==':'; tst(j)='C'; end
         end
         SOLIDNAMES(i,:)=tst;
     end
     
     [C, SPECIATIONNAMES,MASSERR,X]=NRlogX(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5);
end

if flag1==3
     [C, SPECIATIONNAMES,MASSERR,X]=PHREEQCsolve(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,database,acid);
end

end

%---------------__SOLVERS ------NRX NRlogX PHREEQC -----------------------------------

function [C, SPECIATIONNAMES,MASSERR,X]=NRX(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2)

FACTOR=1; LOOP=1;  Xguess=T;%./(1.1^FACTOR)
%loop here to get a decent initial guess based on totals
while LOOP==1
   logC=(KSOLUTION)+ASOLUTION*log10(Xguess);
   C=10.^(logC); % calc species
   R=ASOLUTION'*C-T; show=R./T;
   [value,index]=max(show);
   Xguess(index)=T(index)./(1.1^FACTOR);
   if max(abs(R./T))<1; LOOP=0; end
   FACTOR=FACTOR+0.01;
   if FACTOR>1000; LOOP=0; end %prevent infinite loop
   %pause
end
% %pause

    %Xguess=T./10; 
    TYPX=Xguess;
    [Xguess,masserr,J,RSI,C] = nl_massbalancerrnosolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX);
%     masserr
% pause
    if masserr>1e-4
        Xguess=Xguess/1000000;
        [Xguess,masserr,J,RSI,C] = nl_massbalancerrnosolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX);
    end
%pause
    if masserr>1e-4
        Xguess=Xguess/1000000000000;
        [Xguess,masserr,J,RSI,C] = nl_massbalancerrnosolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX);
    end

    if masserr>1e-4
        Xguess=Xguess/1000000000000000000000;
        [Xguess,masserr,J,RSI,C] = nl_massbalancerrnosolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX);
    end
%masserr
 %   pause

%     X=Xguess; % in case there are no solids. this is the answer.

    %if max(RSI)>0 % only serach for solids if there are some supersat.
        RSI(RSI>0)=min(T); Xguess=[Xguess; RSI]; TYPX=Xguess;
        [X,masserr,J,RSI,C] = nl_massbalancerrsolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2);
    %end

    M=size(X,1); N=size(RSI,1); C=[C; X(M-N+1:M)]; C(C<0)=0;

    SPECIATIONNAMES=[SOLUTIONNAMES
                     SOLIDNAMES];
                 
    %MASSERR=max(abs((100*(masserr./T))));
    
    MASSERR=(abs((100*(masserr./T))));

end

function [C, SPECIATIONNAMES,MASSERR,X]=NRlogX(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,flag3,flag4,flag5)
if flag5==0
% reduce tableau to just the components.

N=length(T);
ASOLUTIONi=ASOLUTION(1:N,:);
KSOLUTIONi=KSOLUTION(1:N);
SOLUTIONNAMESi=SOLUTIONNAMES(1:N);

Xguess=T; TYPX=Xguess; TYPX=ones(size(Xguess));
[Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrnosolid_NR(Xguess,ASOLUTIONi,KSOLUTIONi,ASOLID,KSOLID,T,TYPX,flag2,flag3);
    if flag3==1; show=abs(max(masserr))
    end

% add a row at a time and see if the error is still small

if flag4==1
 for i=N:size(ASOLUTION,1)
     ASOLUTIONi=ASOLUTION(1:i,:);
     KSOLUTIONi=KSOLUTION(1:i);
     SOLUTIONNAMESi=SOLUTIONNAMES(1:i);
     
     %TYPX=Xguess; 
     [Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrnosolid_NR(Xguess,ASOLUTIONi,KSOLUTIONi,ASOLID,KSOLID,T,TYPX,flag2,flag3);
     if flag3==1; show=abs(max(masserr))
     end
 end
end

[Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrnosolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2,flag3);
Xsolution=Xguess; M=length(Xsolution); %N=length(SI);

save Xsolutionguess.mat Xguess

RSIsaved=RSI;

RSIsaved(RSIsaved>0)=0; RSI=RSIsaved;
XguessN=[(Xguess); RSI];
TYPX=[log10(Xguess); ones(size(RSI))]; TYPX=ones(size(TYPX));
end

if flag5==1; load Xsolidsguess.mat; XguessN=X; TYPX=ones(size(X)); end

[X,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(XguessN,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2,flag3);
if flag3==1; max(abs(masserr))    
end
%pause

if max((abs(masserr)))>1e-7 % change pH to make a better initial guess
    %TYPX=T./10;
    pH=-1*KSOLUTION(1); pe=-1*KSOLUTION(2); pH=pH-0.2
    load originaltableau.mat
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pH(oKSOLUTION,oKSOLID,oASOLUTION,oASOLID,pH);
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pe(rKSOLUTION,rKSOLID,rASOLUTION,rASOLID,pe);
    [Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(XguessN,rASOLUTION,rKSOLUTION,rASOLID,rKSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
    [X,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
end

if max((abs(masserr)))>1e-7 % change pH to make a better initial guess
    %TYPX=T./10;
    pH=-1*KSOLUTION(1); pe=-1*KSOLUTION(2); pH=pH-0.4
    load originaltableau.mat
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pH(oKSOLUTION,oKSOLID,oASOLUTION,oASOLID,pH);
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pe(rKSOLUTION,rKSOLID,rASOLUTION,rASOLID,pe);
    [Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(XguessN,rASOLUTION,rKSOLUTION,rASOLID,rKSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
    [X,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
end

if max((abs(masserr)))>1e-7 % change pH to make a better initial guess
    %TYPX=T./10
    pH=-1*KSOLUTION(1); pe=-1*KSOLUTION(2); pH=pH-0.6
    load originaltableau.mat
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pH(oKSOLUTION,oKSOLID,oASOLUTION,oASOLID,pH);
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pe(rKSOLUTION,rKSOLID,rASOLUTION,rASOLID,pe);
    [Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(XguessN,rASOLUTION,rKSOLUTION,rASOLID,rKSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
    [X,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
end

if max((abs(masserr)))>1e-7 % change pH to make a better initial guess
    %TYPX=T./10;
    pH=-1*KSOLUTION(1); pe=-1*KSOLUTION(2); pH=pH-0.8
    load originaltableau.mat
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pH(oKSOLUTION,oKSOLID,oASOLUTION,oASOLID,pH);
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pe(rKSOLUTION,rKSOLID,rASOLUTION,rASOLID,pe);
    [Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(XguessN,rASOLUTION,rKSOLUTION,rASOLID,rKSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
    [X,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
end

if max((abs(masserr)))>1e-7 % change pH to make a better initial guess
    %TYPX=T./10;
    pH=-1*KSOLUTION(1); pe=-1*KSOLUTION(2); pH=pH-1
    load originaltableau.mat
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pH(oKSOLUTION,oKSOLID,oASOLUTION,oASOLID,pH);
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pe(rKSOLUTION,rKSOLID,rASOLUTION,rASOLID,pe);
    [Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(XguessN,rASOLUTION,rKSOLUTION,rASOLID,rKSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
    [X,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
end

if max((abs(masserr)))>1e-7 % change pH to make a better initial guess
    %TYPX=T./10;
    pH=-1*KSOLUTION(1); pe=-1*KSOLUTION(2); pH=pH-1.5
    load originaltableau.mat
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pH(oKSOLUTION,oKSOLID,oASOLUTION,oASOLID,pH);
    [rKSOLUTION,rKSOLID,rASOLUTION,rASOLID]=get_equilib_fixed_pe(rKSOLUTION,rKSOLID,rASOLUTION,rASOLID,pe);
    [Xguess,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(XguessN,rASOLUTION,rKSOLUTION,rASOLID,rKSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
    [X,masserr,J,RSI,C] = NRlogXnl_massbalancerrsolid_NR(Xguess,ASOLUTION,KSOLUTION,ASOLID,KSOLID,T,TYPX,flag2,flag3);
    max(abs(masserr))
    %pause
end

save Xsolidsguess.mat X

M=length(X); N=length(RSI); C=[C; X(M-N+1:M)]; C(C<0)=0;

SPECIATIONNAMES=[SOLUTIONNAMES
                     SOLIDNAMES];
                 
%MASSERR=max(abs((100*(masserr./T))));
%MASSERR=(abs((100*(masserr./T)))); % get rid of relative error.  buggers up on small concs
MASSERR=abs(masserr);

end

function [C, SPECIATIONNAMES,MASSERR,X]=PHREEQCsolve(KSOLID,ASOLID,SOLIDNAMES,KSOLUTION,ASOLUTION,SOLUTIONNAMES,T,flag1,flag2,database,acid);

%---------------------------------------------------------------------------------
%database=['llnl_nosolubleAgCl.dat'];
%minerals=[{'Chlorargyrite'}]
%minerals=num2cell(SOLIDNAMES,[1 2]);
for i=1:size(SOLIDNAMES,1)
    name=SOLIDNAMES(i,:);  name(name == ' ') = [];
    minerals(i)=num2cell(name,[1 2]);
end
minerals=minerals';

%totalnames=[{'Ag '}; { 'Cl '}; {'Na'}; {'N(5)'}]
%totalvector=[AgT    ClT  NaT NO3T];
totalvector=T';
%speciesexport=[{'Ag+'}; {'Cl-'};]
c=0;
for i=3:size(SOLUTIONNAMES,1)
    c=c+1; name=SOLUTIONNAMES(i,:);
    name(name == ' ') = [];
    speciesexport(c)=num2cell(name,[1 2]);
end
speciesexport=speciesexport';

pH=-1*KSOLUTION(1); pe=-1*KSOLUTION(2);
[N,M]=size(ASOLUTION); c=0; 
for i=3:2+M
    c=c+1;
    tst=SOLUTIONNAMES(i,:);
%     for j=1:length(tst)
%         tester=isstrprop(tst(j), 'alpha');
%         if tester==0; tst(j)=' '; end
%     end

   species=['PO4-3']; n=length(species);
   if tst(1:n)=='PO4-3'
       tst1=['P']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['Na']; total2vector(c)=T(c)*2; total2names(c)=num2cell(tst2,[1 2]);
   end

   species=['Ca+2']; n=length(species);
   if tst(1:n)=='Ca+2'
       tst1=['Ca']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['Na']; total2vector(c)=T(c)*2; total2names(c)=num2cell(tst2,[1 2]);
   end
   
  species=['CO3-2']; n=length(species);
   if tst(1:n)=='CO3-2'
       tst1=['C(4)']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['Na']; total2vector(c)=T(c)*2; total2names(c)=num2cell(tst2,[1 2]);
   end
   
   species=['Fe+3']; n=length(species);
   if tst(1:n)=='Fe+3'
       tst1=['Fe'];  totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['Cl']; total2vector(c)=T(c)*3;  total2names(c)=num2cell(tst2,[1 2]);
   end

   species=['Ag+']; n=length(species); 
   if tst(1:n)=='Ag+'
       tst1=['Ag']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['N(5)']; total2vector(c)=T(c)*1;  total2names(c)=num2cell(tst2,[1 2]);
   end
      
   species=['Na+']; n=length(species); 
   if tst(1:n)=='Na+'
       tst1=['Na']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['N(5)']; total2vector(c)=T(c)*1;  total2names(c)=num2cell(tst2,[1 2]);
   end
   
    species=['NO3-']; n=length(species); 
   if tst(1:n)=='NO3-'
       tst1=['N(5)'];  totalnames(c)=num2cell(tst1,[1 2]);
   end
   
   species=['Cl-']; n=length(species);
   if tst(1:n)=='Cl-'
       tst1=['Cl']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['Na']; total2vector(c)=T(c)*1; total2names(c)=num2cell(tst2,[1 2]);
   end
      
   species=['S(-2)']; n=length(species);
   if tst(1:n)=='S(-2)'
       tst1=['S']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['K']; total2vector(c)=T(c)*2; total2names(c)=num2cell(tst2,[1 2]);
   end

   species=['SO4-2']; n=length(species);
   if tst(1:n)=='SO4-2'
       tst1=['S']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['K']; total2vector(c)=T(c)*2; total2names(c)=num2cell(tst2,[1 2]);
   end


   species=['Eu+3']; n=length(species);
   if tst(1:n)=='Eu+3'
       tst1=['Eu']; totalnames(c)=num2cell(tst1,[1 2]);
       %tst2=['K']; total2vector(c)=T(c)*2; total2names(c)=num2cell(tst2,[1 2]);
   end

end

totalnames=[totalnames'];
    %total2names'];

totalvector=[totalvector];%   total2vector];

temp=25;

%-----------------------------------------------------------------------

[solutionspeciesconcs, solutionspeciesnames, solidconcs, solidnames]= ...
callPHREEQC(totalnames,totalvector,pH,pe,temp,minerals,speciesexport,database,acid,flag2);

%parse the solution into variable for each species requested

for j=1:size(solutionspeciesconcs,1)
            name=cell2mat(solutionspeciesnames(j));
            %clean up name so it can be a matlab variable
%             opens = name == '(m';
%             closes = name == 'w)';
%             nestingcount = cumsum(opens - [0 closes(1:end-1)]);
%             name = name(nestingcount == 0);
            name=regexprep(name,'(mol/kgw)','');
            name=regexprep(name,'(eq/kgw)','');
            name=regexprep(name,'[m_]','');
            name = strrep(name,'+','plus');
            name = strrep(name,'-','minus');
            name=regexprep(name,'(','');
            name=regexprep(name,')','');
            %assign variable to corresponding concentration
            txt=[name,'=solutionspeciesconcs(j);'];
            eval(txt)
            %pause
end

indx=ones(size(solidconcs));
for i=1:size(solidconcs,1)
     name=cell2mat(solidnames(i));
     if name(1:2)=='d_' ; indx(i)=0; end
end
c=0;
for i=1:size(solidconcs,1)
    name=cell2mat(solidnames(i));
    if indx(i)==1
        c=c+1;
        SOLIDconcs(c)=solidconcs(i);
        SOLIDnames(c)=solidnames(i);
    end
end


% for j=1:size(SOLIDconcs,1)
%             name=cell2mat(solidnames(j));
%             %clean up name so it can be a matlab variable
%             if name(1:2)=='d_' 
%                 name=regexprep(name,'[d_]','');
%                  name=regexprep(name,':','');
%                 name=regexprep(name,'(','');
%                 name=regexprep(name,')','');
%                 %assign variable to corresponding concentration
%                 txt=[name,'=solidconcs(j);']; eval(txt);
%             end
% end

% outputs

C=solutionspeciesconcs; C(1)=10^-C(1); C(2)=10^(-C(2));

C=[C
    SOLIDconcs'];

X=C(3:2+M);

SUPERA=[ASOLUTION
    ASOLID];
[N,M]=size(SUPERA);

Tcalc=zeros(1,M);
for i=1:M
    for j=1:N
    Tcalc(i)=Tcalc(i)+C(j)*SUPERA(j,i);
    end
end
% tst=[T Tcalc']
% C
% pause
relERR=abs(100*((T-Tcalc')./T));
%MASSERR=max(relERR);
MASSERR=(relERR);

for i=1:size(SOLUTIONNAMES,1)
         tst=SOLUTIONNAMES(i,:);
         for j=1:length(tst)
             if tst(j)=='+'; tst(j)='p'; end
             if tst(j)=='-'; tst(j)='m'; end
             if tst(j)=='('; tst(j)='L'; end
             if tst(j)==')'; tst(j)='R'; end
             if tst(j)==':'; tst(j)='C'; end
         end
         SOLUTIONNAMES(i,:)=tst;
     end
     
     for i=1:size(SOLIDNAMES,1)
         tst=SOLIDNAMES(i,:);
         for j=1:length(tst)
             if tst(j)=='('; tst(j)='L'; end
             if tst(j)==')'; tst(j)='R'; end
             if tst(j)==':'; tst(j)='C'; end
         end
         SOLIDNAMES(i,:)=tst;
     end
     
     %SOLIDNAMES
SPECIATIONNAMES=[SOLUTIONNAMES
    (SOLIDNAMES)];
end