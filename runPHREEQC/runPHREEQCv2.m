function [solutionspeciesconcs, speciesnames, SOLIDconcs, SOLIDnames]=...
    runPHREEQCv2(T,pH,pe,HAconc,totalnames,totalvector,minerals,SURFACECONCS,...
    SURFACENAMES,speciesexport,ionicstrength,dummysurfacearea,database,show,acid);

capacitance=(96500)^2/((2*8.314*298)*(-196*log10(ionicstrength)*dummysurfacearea));

NOOFSOLIDS=size(minerals,1);

% Construct the text file to run PHREEQC from MATLAB-------------------
fileID=fopen('runphreeqc.txt','w');
fclose(fileID);
fileID=fopen('runphreeqc.txt','a');

% define solution  -------------------------------------------------------
fprintf(fileID,'SOLUTION 1\n');
fprintf(fileID,['       pe      ' num2str(pe), '\n']); % the redox value
fprintf(fileID,['       pH      ' num2str(pH), '\n']); % the pH condition
fprintf(fileID,['       temp      ' num2str(T), '\n']); % the temperature
fprintf(fileID,'-units mol/kgw\n'); % the unit of the input; usually mol/L is used
% put in the totals
for i=1:size(totalnames,1)
    tf=strcmp('Cl',cell2mat(totalnames(i)));
    if tf==1
        totaltxt=[cell2mat(totalnames(i)),' ', num2str(totalvector(i)), ' charge\n'];
        fprintf(fileID,totaltxt);
    end
    if tf==0
   totaltxt=[cell2mat(totalnames(i)),' ', num2str(totalvector(i)), ' \n'];
   fprintf(fileID,totaltxt);
    end
end
fprintf(fileID,'\n');

% define equilibrium phases (solids, fix pH, pe) -------------------------
fprintf(fileID,'EQUILIBRIUM_PHASES 1\n'); 
 for i=1:NOOFSOLIDS
    phasestxt=[cell2mat(minerals(i)), '   0.0   0\n'];
    fprintf(fileID,phasestxt);
 end
pHfixline=['       Fix_H+ -',num2str(pH),'          ',acid,' 10.0\n']; 
fprintf(fileID,pHfixline);
fprintf(fileID,'-force_equality true\n');
pefixline=['       Fix_pe ',num2str(-1*pe),'          O2\n']; 
fprintf(fileID,pefixline);
fprintf(fileID,'-force_equality true\n');
fprintf(fileID,'\n');

% define surface (HA, or HFO) -------------------------
fprintf(fileID,'SURFACE 1\n'); 
surfaceline=[cell2mat(SURFACENAMES(1)),' ',num2str(SURFACECONCS(1)),' ',num2str(dummysurfacearea),' ',num2str(HAconc),' \n'];
fprintf(fileID,surfaceline);
for i=2:length(SURFACECONCS)
nextline=[cell2mat(SURFACENAMES(i)),' ',num2str(SURFACECONCS(i)),' \n'];
fprintf(fileID,nextline);
end
fprintf(fileID,'-cd_music \n');
lastline=['-capacitance ',num2str(capacitance),' 1e5 \n'];
fprintf(fileID,lastline);
fprintf(fileID,'\n');

% define outputs of model ------------------------------------------------
fprintf(fileID,'SELECTED_OUTPUT\n');
fprintf(fileID,'-file selected.out\n');
fprintf(fileID,'-selected_out true\n');
fprintf(fileID,'-user_punch true\n');
fprintf(fileID,'-high_precision true\n');
fprintf(fileID,'-reset false\n');
fprintf(fileID,'-simulation false\n');
fprintf(fileID,'-state false\n');
fprintf(fileID,'-distance false\n');
fprintf(fileID,'-time false\n');
fprintf(fileID,'-step false\n');
fprintf(fileID,'-ph false\n');
fprintf(fileID,'-pe false\n');
fprintf(fileID,'-reaction false\n');
fprintf(fileID,'-temperature false\n');
fprintf(fileID,'-alkalinity false\n');
fprintf(fileID,'-ionic_strength false\n');
fprintf(fileID,'-water false\n');
fprintf(fileID,'-charge_balance false\n');
fprintf(fileID,'-percent_error false\n');

%fprintf(fileID,'KNOBS\n');
%fprintf(fileID,'-iterations 1500\n');

% species to export
speciestxt=['-molalities '];
for i=1:size(speciesexport,1)
speciestxt=[speciestxt, cell2mat(speciesexport(i)), ' '];
end
speciestxt=[speciestxt, ' \n'];
fprintf(fileID,speciestxt);
equilibriumphases=['-equilibrium_phases ']; %fprintf(fileID,equilibriumphases)
for i=1:NOOFSOLIDS
     equilibriumphases=[equilibriumphases, cell2mat(minerals(i)), '  '];
end
fprintf(fileID,equilibriumphases);
fclose(fileID);
% 
% % run the model -----------------------------------------------------------
str=['system("$HOME/runPHREEQC/phreeqc runphreeqc.txt out.txt ', database,'")'];
% %system('export LD_PRELOAD=/usr/lib/libstdc++.so.6')
% %system('phreeqc runphreeqc.txt out.txt llnl.dat')
% %!phreeqc runphreeqc.txt out.txt llnl.dat
if show==1
eval(str); % output to the screen
end
if show==0
 evalc(str); % so no screen output
end 
% % import the data from the run and prepare to export----------------
% DIFFERFENT IN OCTAVE THAN MATLAB
fid = fopen('selected.out','rt');
hdr = strtrim(regexp(fgetl(fid),'\t','split')); hdr=hdr(1:2*NOOFSOLIDS+length(speciesexport))';
mat = cell2mat(textscan(fid,repmat('%f',1,numel(hdr))));
fclose(fid);
 
out_PHREEQC=mat';
[n,m]=size(out_PHREEQC); out_PHREEQC=out_PHREEQC(1:n-1,:); 
solutionspeciesconcs=out_PHREEQC(1:length(speciesexport),2);
%solutionspeciesconcs=selectedconcs(1:n-2*NOOFSOLIDS);
solidconcs=out_PHREEQC(length(speciesexport)+1:end,2);
speciesnames=(hdr(1:n-2*NOOFSOLIDS,1));
solidnames=(hdr(n-2*NOOFSOLIDS+1:n,1));
 

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


end
