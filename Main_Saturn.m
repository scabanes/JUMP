clear all
close all
%**********************************************************
%          PDS: PARAMETRI DA SCEGLIERE
%**********************************************************
%---------------------------------------------
%-------------- Nomenclature
Name = 'Polar winds Arrate'; %
roots = '/home/simon/Bureau/Esperimento-DICEA/Observations-Cassini-voyager/'; % Root path..
%**********************************************************

% #########################################################################
%                                                        DATI INFORMAZIONI:
% #########################################################################
run([roots,Name,'/InfosFile_Saturn.m'])
load('besselzeros2_C.mat'); 
%---------------------------------------------
%-------------- Frames & Tocm
Tmax=2;%2401; % Time max wanted, nFrames being its maximum value
Itime=1;
nTime = Tmax-Itime+1;
%
Tocm = 1.; % Fattore necessario per avere cm -> 1. if it is already in cm

%---------------------------------------------
%-------------- Grid
Create_Grid_pol_Saturn


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
% 
% SpectralAnalysis_FB
% 
% SpectralAnalysis_FB_Plots


%--------------------------------------------------------------------------
% ->>> PV Monotonization
%--------------------------------------------------------------------------
%% few points here
% 1/ Note that DerivateVorticita.m is not called here but the vorticity comes
%    from Data2RegularGrid.m in the folder: cd([roots,Name])
% 2/ Attention aussi dans la routine il y a une ligne Vz = Pvpolrm.*0; qui
%    permet de ne pas charger la vitesse zonale.
% 3/ In the north pole the sorting algorithm needs to be 'ascend' instead
%    of descend.

% SortingDirection = 'ascend';
% %
% PV_mono

%--------------------------------------------------------------------------
%                                                              ->>> Efluxes
%--------------------------------------------------------------------------
EFluxes_Spectral
% % % % % 
EFluxes_Spectral_plots
