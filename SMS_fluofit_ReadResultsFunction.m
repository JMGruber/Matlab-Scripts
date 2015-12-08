function [lifetimes,amplitudes] = SMS_fluofit_ReadResultsFunction()
% [lifetimes,amplitudes] = SMS_fluofit_ReadResultsFunction(fname)
% 
% loads a Fluofit results file of a global fit --> extract and export tables of amplitudes and lifetimes
% 
% INPUT
%   fname   -   filename
% OUTPUT
%   lifetimes   -  column vector with lifetimes
%   amplitudes  -  array of amplitudes; size(amplitudes) = [number of lifetimes, number of datasets];
% EXPORTED FILES
%   "filename lifetimes.dat"   -  lifetimes
%   "filename amplitudes.dat"  -  amplitudes
%   "filename summary.dat"     -  all fit parameters
% 
% Note that existing output files will be overwritten without warning
% 
% Bart van Oort, VU University Amsterdam, 16-09-2014

fname='results';
path = uigetdir;
fid = fopen(strcat(path,'\',fname));

i=0;
j1 = 0;
j2 = 0;
j3 = 0;
j4 = 0;
j5 = 0;
IsGlobPar = 0;
IsLocPar = 0;
IsAvgTau = 0;
IsFractInt = 0;
IsFractAmp = 0;
clear GlobParStr LocParStr AvgTauStr FractIntStr FractIntAmp line
while ~feof(fid) % until end of file
    line = fgetl(fid); %# read line by line
    if strcmp(line,'Global Parameters*******************************')
%         disp(line)
        IsGlobPar = 1;
    end
    if strcmp(line,'Local Parameters:*******************************')
%         disp(line)
        IsGlobPar = 0;
        IsLocPar = 1;
    end
    if strcmp(line,'Average Lifetime:')
%         disp(line)
        IsLocPar = 0;
        IsAvgTau = 1;
    end
    if strcmp(line,'Fractional Intensities of the Positive Decay Components:')
%         disp(line)
        IsAvgTau = 0;
        IsFractInt = 1;
    end
    if strcmp(line,'Fractional Amplitudes of the Positive Decay Components:')
%         disp(line)
        IsFractInt = 0;
        IsFractAmp = 1;
    end
    if strcmp(line,'Fitted Decay and Exponential Components:')
%         disp(line)
        IsFractAmp = 0;
    end
    if IsGlobPar==1
        j1 = j1+1;
        GlobParStr{j1} = line;
    end
    if IsLocPar==1
        j2 = j2+1;
        LocParStr{j2} = line;
    end
    if IsAvgTau==1
        j3 = j3+1;
        AvgTauStr{j3} = line;
    end
    if IsFractInt==1
        j4 = j4+1;
        FractIntStr{j4} = line;
    end
    if IsFractAmp==1
        j5 = j5+1;
        FractIntAmp{j5} = line;
    end
end

fclose(fid);

% display results
for i=1:length(GlobParStr); disp(GlobParStr{i}); end
for i=1:length(LocParStr); disp(LocParStr{i}); end
for i=1:length(AvgTauStr); disp(AvgTauStr{i}); end
for i=1:length(FractIntStr); disp(FractIntStr{i}); end
for i=1:length(FractIntAmp); disp(FractIntAmp{i}); end

%% make table of lifetimes
nrLifetimes = strfind(GlobParStr,'t');
nrLifetimes = numel([nrLifetimes{:}])-3 ; % nr. of lifetimes
lifetimes = zeros(nrLifetimes,1);
j=0;
for i=1:length(GlobParStr)
    tmp = strfind(GlobParStr{i},'='); % get pointer to "="
    if ~isempty(tmp) && isempty(strfind(GlobParStr{i},'IRF')) % if a line contains "=", but not "IRF"
        j=j+1;
        lifetimes(j) = str2num(GlobParStr{i}(tmp+1:end-3));
    end
end
lifetimes = flipud(lifetimes);
% disp(lifetimes)

%% make table of amplitudes
nrDataSets = strfind(LocParStr,'Data');
nrDataSets = numel([nrDataSets{:}]);  
amplitudes = zeros(nrLifetimes,nrDataSets);
j=0;
k=0;
for i=1:length(LocParStr)
%     disp(LocParStr{i})
    % if start of new dataset
    if ~isempty(strfind(LocParStr{i},'Data'))
        j=j+1;
        k=0;
    end
    tmp = strfind(LocParStr{i},'='); % get pointer to "="
    if ~isempty(tmp) && isempty(strfind(LocParStr{i},'Bkgr')) && isempty(strfind(LocParStr{i},'Scat')) % if a line contains "=", but not "Bckgr" or "Scat"
        k=k+1;
        amplitudes(k,j) = str2num(LocParStr{i}(tmp+1:end-5));
    end
end
amplitudes = flipud(amplitudes);
% disp(amplitudes)

%% export summary file
fnameout = [fname(1:end-4) ' - summary.dat'];
fidout = fopen(fnameout,'w');
formatSpec = '%s\n';

for i = 1:length(GlobParStr)
    fprintf(fidout,formatSpec,GlobParStr{i});
end
for i = 1:length(LocParStr)
    fprintf(fidout,formatSpec,LocParStr{i});
end
for i = 1:length(AvgTauStr)
    fprintf(fidout,formatSpec,AvgTauStr{i});
end
for i = 1:length(FractIntStr)
    fprintf(fidout,formatSpec,FractIntStr{i});
end
for i = 1:length(FractIntAmp)
    fprintf(fidout,formatSpec,FractIntAmp{i});
end
fclose(fidout);

%% export tables of lifetimes and amplitudes
fnameout = [fname(1:end-4) ' - lifetimes.dat'];
% if exist(fnameout); errordlg(['export file already exists :     ' fnameout],'file'); end
dlmwrite(fnameout,lifetimes,'\t')
fnameout = [fname(1:end-4) ' - amplitudes.dat'];
% if exist(fnameout); errordlg(['export file already exists :     ' fnameout],'file'); end
dlmwrite(fnameout,amplitudes,'\t')
