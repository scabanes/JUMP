clear all
close all
%**********************************************************
%          PDS: PARAMETRI DA SCEGLIERE
%**********************************************************
%---------------------------------------------
%-------------- Nomenclature
Name = 'Full-2D-data/Statistical-JUMP'; %
roots = '/home/simon/Bureau/Esperimento-DICEA/Observations-Cassini-voyager/Jupiter/';
%**********************************************************

% #########################################################################
%                                                        DATI INFORMAZIONI:
% #########################################################################
run([roots,Name,'/InfosFile_JupObs.m'])
load('besselzeros2_C.mat'); 
% load([roots,Name,'/infos.mat'])
%---------------------------------------------
%-------------- Frames & Tocm
nFrames = 3%((endi-starti)/step + 1);
Tmax=nFrames;
Itime=1;
nTime = Tmax-Itime+1;
%
% Tocm = 100.; % Fattore necessario per convertire in cm

%---------------------------------------------
%-------------- Grid
% Create_Grid_pol_Torino_cm

% u = ncread([roots,Name],'u');
% v = ncread([roots,Name],'v');
% imagesc(u(:,:,1,1)')

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
% % EFluxes_Spectral
% % % % 
% % EFluxes_Spectral_plots
