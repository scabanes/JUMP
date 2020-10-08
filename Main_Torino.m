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
% % % % load([roots,Name,'/infos.mat'])
%---------------------------------------------
%-------------- Frames & Tocm

nFrames = 755 %((endi-starti)/step + 1);
Tstep = 1;
Tmax=nFrames;
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
% 
%--------------------------------------------------------------------------
% Waring!! The grid vector r and theta are modified when running the
% programm SpectralAnalysis_FB. Spectral Analysis and PV Monotonization
% have then to be run separately.
%--------------------------------------------------------------------------
%                               ->>> Spectral Analysis in Polar coordinate
%--------------------------------------------------------------------------
% SpectralAnalysis_FB
% 
% SpectralAnalysis_Plots

%--------------------------------------------------------------------------
%                           ->>> Spectral Analysis in Cartesian coordinate
%--------------------------------------------------------------------------
% Here the the velocity field are U(Nx,Ny).
% FieldFrom_Torino as to be choosen in SpectralAnalysis_FF.m 
% (lines 22 & 45-50)
SpectralAnalysis_FF
%   
SpectralAnalysis_FF_Plots

%--------------------------------------------------------------------------
%                                                   ->>> PV Monotonization
%--------------------------------------------------------------------------
%DerivateVorticita
% PV_mono

%--------------------------------------------------------------------------
%                                                              ->>> Efluxes
%--------------------------------------------------------------------------
% EFluxes_Spectral_Torino
% % 
% EFluxes_Spectral_Torino_plots
