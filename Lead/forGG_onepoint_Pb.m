clear; format short e


flag1=2; %flag1 for concentration based search vector (1) or for log based search vector (2).
flag2=1; %flag2 for numerical derivatives (2) or analtyical derivatives (1)
database=[]; 
flag3=0; % flag3 1 for show warnings. 0 no warnings
flag4=0; %if=1 solve tableau one line at a time to get a good initial guess (no solids). really challenging probs only
flag5=0; %if=1 use stored initial guess for including solids

PbT=[1e-7 1e-6 1e-6 1e-7 1e-7 1e-7 1e-7 1e-7 1e-6]; ClT=[0 0 0.1 0.1 0.025 0.025 0 0 0.2]; 
CT=[0.1 0 0 0 0 0 0 0 0]; AcT=[0 0.001 0 0 0 0 0.1 0.1 0]; 
YT=[0 0 0 1e-8 1e-8 1e-7 1e-8 1e-7 1e-7];
pH=[6.5 6.5 6.5 6.5 6.5 7.5 7.5 7.5 7.5]; pe=20.75-pH; 

for i=1:length(PbT)
    
    % make sure totals are in the same order as the Tableau!
    inorganicTOTALS=[PbT(i) ClT(i) CT(i) AcT(i) YT(i)];
    inorganicTOTALS(inorganicTOTALS==0)=1e-16; % get rid of zero values, b/c div by zero error
    TOTALS=[inorganicTOTALS]; %TOTALS=inorganicTOTALS;

    %run the model
    IS=max(inorganicTOTALS);
    if IS==0.001; [Pb(i),PbOH(i),PbCl(i),Pbsolids(i),MASSERR(i)]=Pbtableau001(pH(i),pe(i),TOTALS',flag1,flag2,flag3,flag4,flag5,database); end
    if IS==0.025; [Pb(i),PbOH(i),PbCl(i),Pbsolids(i),MASSERR(i)]=Pbtableau025(pH(i),pe(i),TOTALS',flag1,flag2,flag3,flag4,flag5,database); end
    if IS==0.1; [Pb(i),PbOH(i),PbCl(i),Pbsolids(i),MASSERR(i)]=Pbtableau1(pH(i),pe(i),TOTALS',flag1,flag2,flag3,flag4,flag5,database); end
    if IS==0.2; [Pb(i),PbOH(i),PbCl(i),Pbsolids(i),MASSERR(i)]=Pbtableau2(pH(i),pe(i),TOTALS',flag1,flag2,flag3,flag4,flag5,database); end
   
end

GGPbvalues=[9.23e-9 7.12e-7 3.91e-7 3.52e-8 5.75e-8 3.32e-12 2.37e-9 1.17e-12 2.19e-7];
GGPbOHvalues=[3.7e-10 5.12e-8 1.51e-8 1.36e-9 2.93e-9 NaN 9.16e-10 NaN 7.38e-8];
GGPbClvalues=[NaN NaN 5.2e-7 4.68e-8 2.84e-8 NaN NaN NaN 6.0447e-7];
figure(1); clf
plot(Pb,GGPbvalues,'ko')
hold on
plot(PbOH,GGPbOHvalues,'k*')
plot(PbCl,GGPbClvalues,'ks')
plot(Pb,Pb)

Pbsolids

figure(2); plot(MASSERR)

% calculate kobs
mPb=9.09e4; mPbOH=1.59e6; mPbCl=1.03e6;
kobscalc=mPb*Pb+mPbOH*PbOH+mPbCl*PbCl

kobsmeas1=[7.93e-4 3 2.44e-5 2.29e-3 2.91 7.87e-2 NaN 2.11e-3 4.10e-2];
kobsmeas2=[8.71e-4 2.86 NaN NaN NaN NaN 5.17e-3 5.99e-3 NaN];
kobsmeas3=[6.87e-4 NaN 4.08e-2 1.38e-2 2.70e-1 5.35e-2 6e-3 4.42e-3 4.57e-2];

kobsmeas4=[8.7440e-04   5.0000e-04   5.0000e-04   1.7473e-02   2.1183e-03   9.1191e-02   6.6546e-03   6.5518e-03   5.0000e-04];

figure(3); h=loglog(Pb,kobsmeas4,'ko'); 
set(gca,'linewidth',2,'fontsize',12)
set(h,'markersize',9,'markerfacecolor','b')
xlabel('kobs meas'); ylabel('kobs calc')

[index,value]=sort(kobscalc)