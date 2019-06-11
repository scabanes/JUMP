clear all
close all
%**********************************************************
%          PDS: PARAMETRI DA SCEGLIERE
%**********************************************************
%---------------------------------------------
%-------------- Nomenclature
Name = 'ESP_31'; %
roots = '/media/simon/simon/'; % Root path..
%**********************************************************

% #########################################################################
%                                                        DATI INFORMAZIONI:
% #########################################################################
run([roots,Name,'/InfosFile_Vasca.m'])
load('besselzeros2_C.mat'); 
%---------------------------------------------
%-------------- Frames & Tocm
Tmax=10; % Time max wanted, nFrames being its maximum value
Itime=1;
nTime = Tmax-Itime+1;
%
Tocm = 1.; % Fattore necessario per avere cm -> 1. if it is in cm

%---------------------------------------------
%-------------- Grid
Create_Grid_pol_Vasca


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
% DerivateVorticita
% 
% PV_mono