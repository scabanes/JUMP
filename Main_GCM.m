clear all
close all
%**********************************************************
%          PDS: PARAMETRI DA SCEGLIERE
%**********************************************************
%---------------------------------------------
%-------------- Nomenclature
Name = 'GCM'; %
roots = '/media/simon/simon/DYNAMICO - ';
% run([roots,Name,'/InfosFile_SaturnGCM.m'])
run([roots,Name,'/InfosFile_SaturnGCM_Polar.m'])
CartGrid_Name = '/Data2Regular_CartesianGrid_infos_PN.mat'; % grid and data may change between north and south pole (PN and PS)
%**********************************************************

% #########################################################################
%                                                        DATI INFORMAZIONI:
% #########################################################################
load('besselzeros2_C.mat'); 
%---------------------------------------------
%-------------- Frames & Tocm
nFrames = 2%((endi-starti)/step + 1);
Tmax=nFrames;
Itime=1;
nTime = Tmax-Itime+1;
%
% Tocm = 100.; % Fattore necessario per convertire in cm
% #########################################################################
%                                                                     GRID:
% #########################################################################
Create_Grid_mtx
% #########################################################################
%                                                                   FIELDS:
% #########################################################################
it=0;
% FieldFrom_GCM
FieldFrom_mtx
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
% (lines 22 & 54-57)
% SpectralAnalysis_FF
% %   
% SpectralAnalysis_FF_Plots

%--------------------------------------------------------------------------
%                                                   ->>> PV Monotonization
%--------------------------------------------------------------------------
%DerivateVorticita
% PV_mono

%--------------------------------------------------------------------------
%                                                              ->>> Efluxes
%--------------------------------------------------------------------------
EFluxes_Spectral
% % % % 
EFluxes_Spectral_plots
