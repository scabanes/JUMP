% MATRIX = LOADMTX(FILENAME)
%
% Carica la matrice MATRIX dal un file MTX
% vedi anche SAVEMTX
% by mm & gq

function [a,N]=loadmtx(filename)
filename=sprintf('%s.mtx',filename);
fid=fopen(filename,'rb');
N=fread(fid,[1 2],'ulong')
a=fread(fid,N,'float32');
fclose(fid);
