function [returnOK] = Maps(Vtheta_tot,Vr_tot,Grid_Xp_cm,Grid_Yp_cm,GridR_2C,GridT_2C,Nti,Nri,Nrmin,Nrmax,Ntmin,Ntmax,theta,time)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Vtheta = reshape(Vtheta_tot(:,time),Nti,Nri);
Vr = reshape(Vr_tot(:,time),Nti,Nri);

theta(end)*180./pi

figure; title(['We work on a ',num2str(theta(end)*180./pi),' degres sector']); hold on
contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta,250,'LineStyle','none')
colorbar

MT = nan(size(Grid_Xp_cm));
MT(Ntmin, Nrmin:Nrmax) = 0;
MT(Ntmax, Nrmin:Nrmax) = 0;

h=pcolor(Grid_Xp_cm,Grid_Yp_cm,MT);
set(h, 'EdgeColor', 'none');
% caxis([-1.5 1.5])
% cut
Vtheta=Vtheta(Ntmin:Ntmax,Nrmin:Nrmax);
Vr=Vr(Ntmin:Ntmax,Nrmin:Nrmax);

% Cartesian grid
Grid_Xp_cm=Grid_Xp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
Grid_Yp_cm=Grid_Yp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
% Locally cartesian grid see Read. 2015
GridR_2C=GridR_2C(Ntmin:Ntmax,Nrmin:Nrmax);
GridT_2C=GridT_2C(Ntmin:Ntmax,Nrmin:Nrmax);

figure; title('section we are working on')
contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta,250,'LineStyle','none')
colorbar

figure; title('section on locally cartesian grid')
contourf(GridT_2C,GridR_2C,flip(Vtheta),250,'LineStyle','none')
colorbar
xlabel('r \theta [m]')
ylabel('r [m]')
% caxis([-1.5 1.5])
returnOK = 'Ok';
end

