close all
clear all
clc

filename = '/data/users/twp11/ITO_MIM/Passive/Input/Tools/Run_number.txt';
fileID = fopen(filename,'w');
tline = fgetl(fileID);
n = str2double(tline);
n = n-1;
fprintf(fileID,n);
fclose(fileID);
