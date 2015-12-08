function [] = export(writedir,parameter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
test=inputname(1);
dlmwrite(strcat(writedir,'\',inputname(2)),parameter,'');

end

