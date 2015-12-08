% This script resolves the intensity levels of intensity-time traces, 
% performs various calculations and gives different figures as output
% 
% Intensities are resolved by using the simple algorithm as described in:
% Krüger, Ilioaia, and Van Grondelle,
% Fluorescence Intermittency from the Main Plant Light-Harvesting Complex: Resolving Shifts between Intensity Levels
% Journal of Physical Chemistry B 115:5071–5082, 2011.
%
% Additions to previous version:
% - option to draw all int levels of all traces
% - option to fix background - often more accurate
% - IntvsTauDensityMap.m -> Draws first Panel as density map
% - IntJumpsDensityPlot.m -> Draws a density map of Int before a jump vs. Int after a jump
% - Calc_TransRates.m; -> Calculates kinetic information, based on a 3-state intensity model
% - IntHist.m -> draws intensity histogram of intensity data points (counts / time resolution)
% - I'll send the file powerlaw.m upon request if you want to do a power-law analysis
%
% Notes: very last int level is not resolved, since a nonblinking event cuts it short (either end of acquisition time or bleaching)
%
% First execute 'matlabpool' to activate parallel computing
% Code written for Parallel Computing: n CPUs will decrease execution time by a factor n.
% (though I have not checked if all permutations work for parallel computing)
%
% Run together with:
% - IntShifts_Algorithm.m
% - IntShifts_binints.m
% - IntShifts_CR.m
% - IntShifts_P4vars.m
% - IntShifts_TestSM.m
% - IntShifts_TrimTrace.m
% - IntvsTauDensityMap.m -> Draws "first Panel" as density map. Change necessary parameters in m file
% - IntJumpsDensityPlot.m -> Draws a density map of Initial Int vs. Final Int. Change necessary parameters in m file.
% - Calc_TransRates.m; -> Calculates kinetic information, based on a 3-state intensity model. Change necessary parameters in m file
% - IntHist2.m -> draws intensity histogram of intensity data points (counts / time resolution). Change necessary parameters in m file. 
% 
% Tjaart Krüger
% Vrije Universiteit Amsterdam
% (c) 2010-2012
%
% Improvements: 
% - option to exclude files
% - add option to average instead of bin and average 1-2-3,2-3-4,3-4-5,...
% - include dark fractions (often identified as photonbursts). photonbursts should be identified on the basis of dwell times + ratio of on:off
% - Don't specify all k's separately but calculate all for 1
% user-defined sigma deviation
% - fprintf: 1...n possible to display?
% - Don't need all temp parameters with variable sizes: calculations cost time.
% - Check testSM and photonburstfactor. First doesn't work well when dim
% fraction is included
% Better to define large vectors/matrices and remove all the zeros in the end!


%% Initialisation
tic;
clear all; close all;
allfiles = 3; % A:B or matrix [] of individual spec files. Use positive integers only!
skipfiles = []; % ditto
%readdir = 'D:\data\Werk\VU\data\Michael\TCSPC';
[filename1, pathname]=uigetfile('*.pt3', 'T3 Mode data:', 0, 0);
writedir = [pathname 'Intanalysis'];
%outfile = [pathname filename(1:end-4) '.out'];
filename2 = 'trace';

ki = [3.25 1.96 1.37 1]; % sigma-deviation for 1..4 consecutive trace points, resp.  Paper: ki=[3.25 1.96 1.37 1] But 10-ms data noisier than "normal" shot noise, prob due to triplets as 1 source
Nf = 1.5; % noise multiplication factor (1.5 - 2 for good S/N)

thrliveI = 3000;  % threshold for final intensity level (in cps, background included, intbin excluded)
thrlivet = 10;  % minimal survival time (in seconds)
thr2c = 15000;  % intensity threshold when definitely >1 complex (in cps, background included, intbin excluded)
intbin = 1;   % number of intensity values to be binned together (input)---- for new version the definition is slightly different
timeres = 0.01; % time resolution (excl intbin)
fixbg = true; % usually more accurate
guessbg = 500; % estimated background; used to test for SM (in cps). If fixbg=false, set rather too large than to small
testSM = false; % test for single quantum unit by <=2-time step into quenched state
testphotonburst = false;
photonburstfactor = 10;   % factor by which dwelling time in Q state exceeds dwelling time in unQ states to define photonburst (typically 10)
intermedfactor = 1; % parameter used to test if there are too many unnatural fluctuations (signifying unstable complex)
                    % decrease if more fluctuations should be incl. (e.g. 0.9)

maxsize = 60/0.01;  % max size of any of the data sets

histbin = 1;  % number of intensity levels to be binned for plots (output)
Qextent = 0.7;  % extent of quenching (for panel 4)

figtitle = '';

drawlevels = true;  % draw intensity levels
drawdensitymaps = false; % draw intensity maps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Access beyond this point at own risk                    %%%                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate matrix sizes
nfiles=length(allfiles);
sortedtraces=zeros(nfiles,5);  % [1..5] = [allusedtraces alldead alldble allphotonburst allusedtracebg]
allnc=zeros(nfiles,1); % number of intensity levels for each complex
alltauon=zeros(maxsize,nfiles);
alltauoff=zeros(maxsize,nfiles);
allinton=zeros(maxsize,nfiles);
%allintoff=ones(maxsize,maxsize)*10*thr2c; % since int can be 0 after bg subtraction
allintoff=ones(maxsize,nfiles)*10*thr2c; % since int can be 0 after bg subtraction
allint=ones(maxsize,nfiles)*10*thr2c;
alldwelltime=zeros(maxsize,nfiles);
allstarttime=zeros(maxsize,nfiles);
AllDelaytimeShiftIndices=zeros(maxsize,nfiles);

alltlive=0;

%% Core
i = 1;
for tr=allfiles
    complex=tr;
    j=1;
    goodfile=true; %i.e. not to be skipped
    skipfiles=sort(skipfiles);
    while (j<=length(skipfiles))&&(skipfiles(j)<=complex)&&(goodfile)
        if complex==skipfiles(j)
            goodfile=false;
        end
        j=j+1;
    end
    if goodfile
        tempindex = zeros(1,5);
    %    readfile = fullfile(readdir,[filename int2str(tr)]);
    %    trace=dlmread(readfile);
    %    trace=dlmread(strcat(readdir,'\trace (',int2str(i),')'));
        [trace,delaytimes] = read_pt3_v4(timeres,pathname,filename1);

        dlmwrite(fullfile(writedir,['trace' int2str(tr)]),trace,' ');
        
        if intbin>1     % bin intensity values
            trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
        end
        temptimeres=timeres;    %*****************not timeres but switch to real #counts
        trace(:,2)=trace(:,2)*temptimeres;
        guessbgc=guessbg*temptimeres*intbin;
        thrliveIc=thrliveI*temptimeres*intbin;

        tlive=IntShifts_TrimTrace(trace,thrliveIc); % trim trace into surviving time
        if tlive<thrlivet/trace(1,1);
            tempindex(2)=tr; 
        else
           trace=IntShifts_CR(trace,thr2c*temptimeres*intbin,tlive); % remove excessively large intensities
           dble=false; photonburst=false;
           if testSM
               SM=IntShifts_testSM(trace,guessbgc*4/3,tlive,thrliveIc);
               if ~SM
                  dble=true;
                  tempindex(3)=tr;
               end
           end
           if (testphotonburst)&&((length(find(trace(1:tlive,2)<thrliveIc))>photonburstfactor*length(find(trace(1:tlive,2)>thrliveIc)))...
                                ||(sum(trace(1:tlive,2)<intermedfactor*thrliveIc & trace(1:tlive,2)>(guessbgc+thrliveIc)/2)... % intermediate intensities
                                > sum(trace(1:tlive,2)>intermedfactor*thrliveIc | trace(1:tlive,2)<(guessbgc+thrliveIc)/2)))       % are more than Q and unQ
                photonburst=true;
                tempindex(4)=tr; 
           end

            if (~photonburst)&&(~dble)
                [intlevels,inttimes,intstart,SM] = IntShifts_Algorithm(trace,ki,Nf,tlive,thr2c*temptimeres*intbin);
                if ~SM
                    tempindex(3)=tr; 
                    dble=true;
                end

                % binned delaytimes is calculated before bg is subtracted
                for j = 1:length(intstart);
                    DelaytimeShiftIndices(j) = sum(trace(intstart(j):intstart(j)+inttimes(j)-1,2));  %these are the sizes, not really the indices - see next command
                end
                DelaytimeShiftIndices = [1 DelaytimeShiftIndices(1:end-1)+1];

                % estimate background
                if fixbg
                    bg=guessbgc;
                else
                    if tlive<trace(end,1)/trace(1,1);
                        blI=trace(tlive:end,2);
                        avgblI=sum(blI)./length(blI);
                        [bg,bgi] = min([intlevels avgblI]);
                    else
                        [bg,bgi] = min(intlevels);
                    end
                    if (bgi<=length(inttimes))&&(inttimes(bgi)*intbin<=2)   % if short dwell times don't represent bad estimations, exchange "2" with "1"
                        levels2 = intlevels;
                        times2 = inttimes;
                        while (bgi>1)&&(bgi<length(levels2))&&(length(levels2)>3)&&(times2(bgi)<=2)
                            levels2 = [levels2(1:bgi-1) levels2(bgi+1:length(levels2))];
                            times2 = [times2(1:bgi-1) times2(bgi+1:length(times2))];
                            [bg,bgi] = min(levels2);
                        end
                    end
                    bg = min(bg,guessbgc);
                end
                intlevels = intlevels-bg;
    %            intlevels(intlevels<-2)=0; % can remove this command, if needed
                inttimes = inttimes.*trace(1,1);  %is this simply times timeres??
                intstart = intstart.*trace(1,1);  %is this simply times timeres??
                alltlive = alltlive+tlive;

                % calculations for 4th panel of figure
                if ~dble
%                     tempalltauon = zeros(maxsize,1);
%                     tempalltauoff = zeros(maxsize,1);
%                     tempallinton = zeros(maxsize,1);
%                     tempallintoff = zeros(maxsize,1);
%                     tempallint = ones(maxsize,1)*10*thr2c;
%                     tempalldwelltime = zeros(maxsize,1);
%                     tempallstarttime = zeros(maxsize,1);

                    P4vars = IntShifts_P4vars(intlevels,inttimes,thrliveIc-bg,Qextent);
                    if ~isempty(P4vars)
                         l = size(P4vars,1);
%                         tempalltauon(1:l,1) = P4vars(1:end,1);     
%                         alltauon(:,i) = tempalltauon;
%                         tempalltauoff(1:l,1) = P4vars(1:end,2);     
%                         alltauoff(:,i) = tempalltauoff;
%                         tempallinton(1:l,1) = P4vars(1:end,3);     
%                         allinton(:,i) = tempallinton;
%                         tempallintoff(1:l,1) = P4vars(1:end,4);     
%                         allintoff(:,i) = tempallintoff;
                        
                        alltauon(1:l,i) = P4vars(:,1);     
                        alltauoff(1:l,i) = P4vars(:,2);     
                        allinton(1:l,i) = P4vars(:,3);     
                        allintoff(1:l,i) = P4vars(:,4);     
                    end

                     l = length(intlevels);
%                     tempallint(1:l,1) = intlevels;
%                     allint(:,i) = tempallint;
%                     tempalldwelltime(1:l,1) = inttimes;
%                     alldwelltime(:,i) = tempalldwelltime;
%                     tempallstarttime(1:l,1) = intstart(1:l);
%                     allstarttime(:,i) = tempallstarttime;

                    allint(1:l,i) = intlevels;
                    alldwelltime(1:l,i) = inttimes;
                    allstarttime(1:l,i) = intstart;
                    AllDelaytimeShiftIndices(1:l,i) = DelaytimeShiftIndices;
                    dlmwrite(fullfile(writedir,['delaytimes_trace' int2str(tr)]),delaytimes',' ');
                    
                    tempindex(1) = tr; 
                    tempindex(5) = bg;
                    allnc(i) = length(intlevels);
                end
            end
        end
        sortedtraces(i,:) = tempindex;
    end
    fprintf(1,'%6.2f\n',tr);
    i=i+1;
end
% remove zeros from matrices
allusedtraces=sortedtraces(sortedtraces(:,1)>0,1);
alldead=sortedtraces(sortedtraces(:,2)>0,2);
alldble=sortedtraces(sortedtraces(:,3)>0,3);
allphotonburst=sortedtraces(sortedtraces(:,4)>0,4);
allusedtracebg=sortedtraces(sortedtraces(:,5)>0,5);
allint=allint(allint<(10*thr2c));
alldwelltime=alldwelltime(alldwelltime>0);
allstarttime=allstarttime(allstarttime>0);
allnc=allnc(allnc>0);
alltauon=alltauon(alltauon>0);
alltauoff=alltauoff(alltauoff>0);
allinton=allinton(allinton>0);
allintoff=allintoff(allintoff<10*thr2c);

%% Plot output
if ~isempty(allint)
    if drawlevels
        readfile = fullfile(writedir,filename2); % int2str(i)]);
        ri=1;
        for tr=1:length(allnc)
            rf=ri+allnc(tr)-1;
%            trace=dlmread(strcat(readfile,[' (' int2str(allusedtraces(tr)) ')']));
            trace=dlmread(strcat(readfile,int2str(allusedtraces(tr))));
            if intbin>1     % bin intensity values
                trace=IntShifts_binints(trace,intbin); % Alternatively, employ IntShifts_avgints.m
            end
            figure; h=plot(trace(:,1),trace(:,2).*timeres-allusedtracebg(tr),'g');
            hold on;
            for k = ri:rf
                x = [allstarttime(k) allstarttime(k)+alldwelltime(k)];
                y = [allint(k) allint(k)];
                plot(x,y,'LineWidth',4,'Color','k');
            end
            ri=ri+allnc(tr);
            saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(tr)) '.jpg']));
            saveas(h,fullfile(writedir,['trace' int2str(allusedtraces(tr)) '.pdf']));
            if mod(tr,10)==0
                close all;
            end
        end
    end

    % Global overview of switches
    h=figure; subplot(2,2,1); semilogx(alldwelltime,allint,'.');
    xlabel('Dwell time (s)'); ylabel(['Intensity (c/',int2str(timeres*1000*intbin),' ms)']);
    title(figtitle);

%    histbin=histbin*timeres; % Take note!
    maxbin=ceil(max(allint)/histbin)*histbin;
    %maxbin=ceil(thr2c*timeres/histbin)*histbin;
    minbin=ceil(min(allint)/histbin)*histbin;

    markers=minbin:histbin:maxbin;
    m=length(markers);
    bin.times=zeros(1,m);
    bin.n2=zeros(1,m);
    for tr=1:length(allint) 
        p=ceil((allint(tr)-minbin)./histbin);
        if (p<=0)||(isnan(p)), p=1; end
%        if p>maxbin, p=maxbin; end
        bin.times(p)=bin.times(p)+alldwelltime(tr);
        bin.n2(p)=bin.n2(p)+1;
    end
    bin.times=bin.times./(alltlive*timeres)*100;
    bin.n2=bin.n2./(alltlive*timeres)*60;
    subplot(2,2,2); bar(markers,bin.n2); 
    xlabel(['Intensity (c/',int2str(timeres*1000*intbin),' ms)']); ylabel('Access frequency (min^{-1})');
    subplot(2,2,3); bar(markers,bin.times);
    xlabel(['Intensity (c/',int2str(timeres*1000*intbin),' ms)']); ylabel('Total dwell time (%)');

    edges=timeres*(exp(0:log(2):9));
    tausortoff=histc(alltauoff,edges)./(alltlive*timeres)*60;
    if ~isempty(tausortoff)
        subplot(2,2,4); semilogx(edges,tausortoff,'.-'); ylabel('Switching frequency (min^{-1})'); xlabel('Dwell time (s)');
    end
    tausorton=histc(alltauon,edges)./(alltlive*timeres)*60;
    if ~isempty(tausorton)
        hold on;
        subplot(2,2,4); semilogx(edges,tausorton,'o-');
    end
    saveas(h,fullfile(writedir,'analysis.jpg'));
    
    if drawdensitymaps
        IntvsTauDensityMap; % First Panel as density map
        IntJumpsDensityPlot; % Draws a density map of Initial Int vs. Final Int
    end
    
    Calc_TransRates;
    
%    IntHist2;  need to change readdir!!
    
end

%% File output
%dlmwrite(strcat(writedir,'\','parameters'),['bg  ',int2str(bg)],'');    
if ~isempty(allusedtraces)
    dlmwrite(strcat(writedir,'\parameters'),['k1 k2 k3 k4 Nf ',int2str(ki(1)),' ',int2str(ki(2)),' ',int2str(ki(3)),' ',int2str(ki(4)),' ',int2str(Nf)],'');
    dlmwrite(strcat(writedir,'\parameters'),['thrliveI thrlivet  ',int2str(thrliveI),' ',int2str(thrlivet)],'-append','delimiter','');
    dlmwrite(strcat(writedir,'\parameters'),['thr2c maxbg ',int2str(thr2c),' ',int2str(guessbg)],'-append','delimiter','');
    dlmwrite(strcat(writedir,'\parameters'),['intbin histbin Qextent  ',int2str(intbin),' ',int2str(histbin),' ',int2str(Qextent)],'-append','delimiter','');

    dlmwrite(strcat(writedir,'\panel1'),[alldwelltime allint],' ');  
    dlmwrite(strcat(writedir,'\panel23'),[markers;bin.n2;bin.times]',' ');
    dlmwrite(strcat(writedir,'\panel4'),[edges' tausortoff tausorton],' ');
    
    dlmwrite(strcat(writedir,'\AllDelaytimeShiftIndices'),AllDelaytimeShiftIndices',' ');

    dlmwrite(strcat(writedir,'\complexes used'),allusedtraces',''); 
    dlmwrite(strcat(writedir,'\complexes dead'),alldead',''); 
    dlmwrite(strcat(writedir,'\complexes double'),alldble',''); 

    save(strcat(writedir,'\analysis.mat'));
end
toc