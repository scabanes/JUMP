clear all
close all
%**********************************************************
%          PDS: PARAMETRI DA SCEGLIERE
%**********************************************************
%---------------------------------------------
%-------------- Nomenclature
%------------------------------------------------------- Polar Forcing
% Name = 'Exp001'; %
% roots = '/media/simon/simon/Polare/'; % Root path..
%------------------------------------------------------- Cartesian Forcing
Name = 'ESP_31_180'; %
roots = '/media/simon/simon/'; % Root path..
%**********************************************************

% #########################################################################
%                                                        DATI INFORMAZIONI:
% #########################################################################
run([roots,Name,'/InfosFile_Vasca.m'])
load('besselzeros2_C.mat'); 
%---------------------------------------------
%-------------- Frames & Tocm
nFrames=3;%1180;%
Tmax=nFrames;%2401; % Time max wanted, nFrames being its maximum value
Itime=1;%1151;
nTime = Tmax-Itime+1;
%
Tocm = 1.; % Fattore necessario per avere cm -> 1. if it is already in cm

%---------------------------------------------
%-------------- Grid
Create_Grid_pol_Vasca % to be replaced by Create_Grid_mtx.m

% #########################################################################
%                                                                   FIELDS:
% #########################################################################
it=0;
%--------------------------------------------------------------------------
% POSSIBLE INPUTS:
% - FieldFrom_2DJupiter_Obs
FieldFrom_Select
% #########################################################################
%                                                                      RUN:
% #########################################################################
%--------------------------------------------------------------------------
% Waring!! The grid vector r and theta are modified when running the
% programm SpectralAnalysis_FB. Spectral Analysis and PV Monotonization
% have then to be run separately.
%--------------------------------------------------------------------------
% ->>> Spectral Analysis
%--------------------------------------------------------------------------

% SpectralAnalysis_FB

% SpectralAnalysis_FB_Plots

%--------------------------------------------------------------------------
%                           ->>> Spectral Analysis in Cartesian coordinate
%--------------------------------------------------------------------------
% SpectralAnalysis_FF
%   
% SpectralAnalysis_FF_Plots


%--------------------------------------------------------------------------
% ->>> PV Monotonization
%--------------------------------------------------------------------------

% DerivateVorticita
% % %
% SortingDirection = 'descend';
% PV_mono

%--------------------------------------------------------------------------
% ->>> Energy & Enstrophy fluxes% % DerivateVorticita
% % % % % % %   
% % PV_mono

%--------------------------------------------------------------------------
%                                                              ->>> Efluxes
%--------------------------------------------------------------------------
EFluxes_Spectral
% % % % 
EFluxes_Spectral_plots


%--------------------------------------------------------------------------
% One has to choose how to discretise the smallest scales:
% dl=0.4;
% EFlux_Energy_old
% EFlux_Enstrophy

% % % passo_fr = 1;
% % % Vl = [1:128];
% % % % 
% % % EFlux
% % % EFlux_Plots