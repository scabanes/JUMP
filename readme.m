% Step by step with JUMP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% VASCA & TORINO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMENTS: One may apply the following chosing Vasca or Torino suffix in 
%% the following codes:
%% Infos:
%->>> InfosFile_suffix.m
% This file has to be nest in the "roots" repository where velocity data of  
% zonal and radial fields are.
% _________________________________________________________________________
%                            A. Spectral analysis
% _________________________________________________________________________
%% MAIN:
%->>> Main_suffix.m
% paths and few infos heve to be entered here and Sub-routines
% appropriately uncommented to run. 
% Note that Create_Grid_pol_suffix.m sets the appropriate grid on which data 
% have been interpolated.
%% Sub-routines:
% 1/ Statistic-----------------------
%->>> SpectralAnalysis_FB.m
% Here we compute the spectral analysis. Everything is in cm !!
% 2/ plots---------------------------
%->>> SpectralAnalysis_Plots.m
% Plots files saved in the Main.
%% Functions:
%->>> FourierBesselDecomp.m
% Function that compute the spectral decomposition.
%->>> Maps.m
% This maps the first time step to see the sector you are working on.
%->>> Windowing.m
% This is a function that creates matrices for signal windowing in a
% choosen direction and with Tukey or Hanning functions.
%->>> loadmtx.m
% It loads .mtx files of the velocity fields.
%->>> lbesselzeros2_C.mat
% Are the required zeros of the Bessel function needed to compute the
% spectral decomposition.
% _________________________________________________________________________
%                    B. Potential Vorticity monotonization
% _________________________________________________________________________
%% MAIN:
%->>> Main_suffix.m
% paths and few infos heve to be entered here and Sub-routines
% appropriately uncommented to run. 
% Note that Create_Grid_pol_suffix.m sets the appropriate grid on which data
% have been interpolated.
%% Sub-routines:
%->>> DerivateVorticita.m
% Leads to relative and potential vorticity for all frames. Data are saved
% in the appropriate folder.
%->>> PV_mono.m:
% PV is monotonized/sorted in order to extract the equivalent Thorpe scale
% denoted L_M. Three averaging procedures can be followed as detailed in
% PV_mono_Torino.
% By looking at the number of points on the radius and in azimuth one might
% have to take out some points. PtDaTogliere = (pt1,pt2,...) give the index
% of the points to take out on the radius when there is too few data to
% have a complete radial profil.
%                           --> L_M is the result of PV mono.
%% Plot-routines:
%->>> PV_Plots
% Plots files saved in the Mains.
% _________________________________________________________________________
%                    C. Energy and Enstrophy Fluxes
% _________________________________________________________________________
% This codes are based on: Evidence for the double cascade scenario in two-
% dimensional turbulence,  G. Boffetta1 and S. Musacchio 2010. 
% DOI: 10.1103/PhysRevE.82.016307
% !! Here data have to be interpolated on a cartesian grid.
%% MAIN:
%->>> Main_suffix.m
% paths and few infos have to be entered here and Sub-routines
% appropriately uncommented to run. 
% Note that Create_Grid_pol_suffix.m sets the appropriate grid on which data
% have been interpolated.
%% Sub-routines:
%->>> EFlux.m:
%-------------------To do:
% Pb de memoire pour accumuler tous les temps. Faire qqu chose.
%--------------------------
% It compute the energy and enstrophy flux using a filtering methode.
%->>> GaussianFilter.m
% The filtering approach is computed using the convolution og the fields
% with a gaussian function G. Field_fil = conv(Field,G). Fluxes are
% obtained by filtering at different scales the velocity and vorticity
% fields.
%% Plot-routines:
%->>> EFlux_Plots.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ## Load data from the Turnbase for Torino experiments:
% https://turbase.cineca.it/init/routes/#/logging/view_dataset/99/tabmeta
% in -> files/data/
