clear all
close all
%**********************************************************
%          PDS: PARAMETRI DA SCEGLIERE
%**********************************************************
%---------------------------------------------
%-------------- Nomenclature
Name = 'EXPT09'; %
roots = '/media/simon/simon/Torino/';
%**********************************************************

% #########################################################################
%                                                        DATI INFORMAZIONI:
% #########################################################################
run([roots,Name,'/InfosFile_Torino.m'])
load('besselzeros2_C.mat'); 
load([roots,Name,'/infos.mat'])
%---------------------------------------------
%-------------- Frames & Tocm
nFrames = ((endi-starti)/step + 1);
Tmax=2;
Itime=1;
nTime = Tmax-Itime+1;
%
Tocm = 100.; % Fattore necessario per convertire in cm

%---------------------------------------------
%-------------- Grid
Create_Grid_pol_Torino_cm


% #########################################################################
%                                                                      RUN:
% #########################################################################
%--------------------------------------------------------------------------
% Waring!! The grid vector r and theta are modified when running the
% programm SpectralAnalysis_FB. Spectral Analysis and PV Monotonization
% have then to be run separately.
%--------------------------------------------------------------------------

% ->>> Spectral Analysis
SpectralAnalysis_FB

%SpectralAnalysis_Plots

% ->>> PV Monotonization
%DerivateVorticita

%PV_mono