function [Vr] = Maps_Galperin2015(Vtheta_tot,Vr_tot,Grid_Xp_cm,Grid_Yp_cm,Nti,Nri,Nrmin,Nrmax,Ntmin,Ntmax,theta,time)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


Vtheta = reshape(Vtheta_tot(:,time),Nti,Nri);
Vr = reshape(Vr_tot(:,time),Nti,Nri);

theta(end)*180./pi

figure; title(['We work on a ',num2str(theta(end)*180./pi),' degres sector']); hold on
contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta,250,'LineStyle','none')

MT = nan(size(Grid_Xp_cm));
MT(Ntmin, Nrmin:Nrmax) = 0.5;
MT(Ntmax, Nrmin:Nrmax) = 0.5;

h=pcolor(Grid_Xp_cm,Grid_Yp_cm,MT);
set(h, 'EdgeColor', 'none');
% cut
Vtheta=Vtheta(Ntmin:Ntmax,Nrmin:Nrmax);
Vr=Vr(Ntmin:Ntmax,Nrmin:Nrmax);

Grid_Xp_cm=Grid_Xp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
Grid_Yp_cm=Grid_Yp_cm(Ntmin:Ntmax,Nrmin:Nrmax);

figure; title('section we are working on')
contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta,250,'LineStyle','none')

end

