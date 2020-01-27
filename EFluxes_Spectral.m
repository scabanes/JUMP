%--------------------------------------------------------------------------
% Alias: test_filtraggio_23Janv20 in dev
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PTS: Parameter To Set in InfosFile_Torino.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is the number of l in pixels on which I calculate the fluxes. It is 
% theoretically limited by Nx and Ny the size of the velocity field in 
% pixels. But it appears that with Nl = Nx & Ny, if the velocity field is
% not periodic the filtering does not manage to filter out every large scales 
% structures. If you impose Nl >> Nx then you finally get to a flux = 0 at 
% largest scales. This might be improved by using a different filter than a
% currently used gaussian, e.g une fonction porte. This can be changed in
% func_Spatial_GaussianFilter.m and/or func_Spectral_GaussianFilter.m
% % % % Nl = 221;
% In spectral filtering it results that l=1 is an uresolved scale that
% corresponds to the lunquist frequency and is theoretically out of
% resolution, consequently the resulting filtered flow is meaning less at l=1. 
% Spatial filtering however manage with it by filtering nothing.
% % % % vl=([2:1:Nl]);
% Here we can c
% % % % FSpatial=0;% set to 1 to have spatial filter and 0 to have spectral filter.
% % % % gaussian = 1; % If gaussian = 1 then G~exp(-6r^2/l^2) if 2 G~exp(-9r^2/(2l^2))
%---------------------------------------------
% Cartesien plein Nx=Ny=150
% idxmin = 40;
% idxmax = 189;
% idymin = 50;
% idymax = 199;
% Cartesien Nx=Ny=221
% % % % idxmin = 1;
% % % % idxmax = 221;
% % % % idymin = 26;
% % % % idymax = 246;
%---------------------------------------------
% Here I used Roland's interpolation through the data merged.civ2. The
% following gapt allow to load data in time.
% % % % gapt = 30;
criterionNaN = 20; % UNSUED -- window on which the mean to replace Nan is calculated
%---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Notes 
%                               &
%                       limit of our approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1/ The fluxes are calcultaed as a function of l set as vl vectore. l=1 is
% non resolved and the smallest scales are poorly resolved. Here l is in 
% pixels. Then fluxes can be presented as a function of k=2pi/l in [k]=px^-1
% Then everything has to be converted in metrics.
% 2/ The biggest issue: When the velocity field is
% non periodic both technics consider strong velocity gradient at the edge 
% of the domain are small scales that have to be filtered out. This
% introduce very strong fluxes when the spatial derivative is computed in
% those outermost region. Consequently, in non periodic field, as it is the
% case for Torino data, we have to average the fluxes in the innermost par
% of the domain and cut out the edge. I see no other solution for that.
Flusso_l = zeros (length([idymin:idymax])*length([idxmin:idxmax]),nFrames,length(vl));
Flusso_Enstro_l = zeros (length([idymin:idymax])*length([idxmin:idxmax]),nFrames,length(vl));
% ##################################################################################################################################################
%                                                                                                                                   INTERPOLAZIONE :
% ##################################################################################################################################################
iit=0;
for t=1:nFrames
    t
% file to load    
file = [roots,Name,'/merged.civ2/fig_',num2str(numt),'.nc'];
% compter les images
numt = numt+gapt;
% Velocity fields of cartesian velocities
UU = double(ncread(file,'U'));
VV = double(ncread(file,'V'));
curl = double(ncread(file,'curl'));
% regular cartesian coordinates
xx = double(ncread(file,'coord_x'));
yy = double(ncread(file,'coord_y'));
%
x = xx(idxmin:idxmax);
y = yy(idymin:idymax);
Nx = length(x);
Ny = length(y);
dx=mean(x(2:end)-x(1:end-1));
dy=mean(y(2:end)-y(1:end-1));
Lx = x(end)-x(1); % en m
Ly = y(end)-y(1); % en m
% Recreat the mesh grid
[XX,YY] = meshgrid(xx,yy);

% Cut to a rectangular table
U = UU(idymin:idymax,idxmin:idxmax);
V = VV(idymin:idymax,idxmin:idxmax);
vort = curl(idymin:idymax,idxmin:idxmax);
X = XX(idymin:idymax,idxmin:idxmax);
Y = YY(idymin:idymax,idxmin:idxmax);
% Fill in the Nans in the table with 0 values
U(isnan(U))=0.;
V(isnan(V))=0.;
vort(isnan(vort))=0.;
% U = fillmissing(U,'movmedian',criterionNaN);  
% V = fillmissing(V,'movmedian',criterionNaN);  
% vort = fillmissing(vort,'movmedian',criterionNaN);
if(0==1)
 figure
 contourf(X,Y,U, 100, 'LineStyle','none');
 colorbar
 size(U)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                 Field spectral decomposition:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          FOURIER DECOMPOSITION for U:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% la fft deve essere normalisata per il numero di punti in azimuth che corrisponde
% al numero dei modi Nm che escono della fft. Tra questi Nm modi ci sono M modi
% zonali m che corrisponde al numero dei punti unici simetrici per M ne l vettore _fft.
%-----------U
U_fft=fft(U,[],2)/Nx; % !! it is normalized
V_fft=fft(V,[],2)/Nx; % !! it is normalized
Nkx=ceil((Nx+1)/2);
%
U_2fft=fft(U_fft,[],1)/Ny; % !! it is normalized
V_2fft=fft(V_fft,[],1)/Ny; % !! it is normalized
Nky=ceil((Ny+1)/2);
% vecteur des modes kx
kx = zeros(1,Nx);
kx(1,1:Nkx) = [1:Nkx]-1;
if(mod(Nx,2)==0) % pair
kx(1,Nkx+1:Nx)=-(Nkx-2:-1:1);
else % impair
kx(1,Nkx+1:Nx)=-(Nkx-1:-1:1);
end
% vecteur des modes ky
ky = zeros(Ny,1);
ky(1:Nky) = [1:Nky]-1;
if(mod(Nx,2)==0) % pair
ky(Nky+1:Ny)=-(Nky-2:-1:1);
else % impair
ky(Nky+1:Ny)=-(Nky-1:-1:1);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Energy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pour l energie on a 0.5 * U^2
% E = 0.5.* (U_2fft.*conj(U_2fft) + V_2fft.*conj(V_2fft));
E = 0.5.*U_2fft.*conj(U_2fft);
E_Tky = sum(E(:,:),2);
%
E_ky = zeros([Nky,1]);
% ce vecteur est juste pour organiser la collecte de modes.
vky = ky.*sign(ky); 
for iky_p=0:Nky-1 % que les ky positifs
    vv=(vky==iky_p);
    E_ky(iky_p+1) = sum(E_Tky(vv)); % j'aurai pu ajouter simplement un *2
end
if(0==1)
figure
loglog(kx(1:Nkx).*2.*pi./Nx,E(10,1:Nkx))
hold on
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Uncomment to use
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This part is to compute some test on Spatial and Spectral filtering
% approaches without going into the associated function. Uncomment to use..
% % % % % % % l=4;
% % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % Spatial_GaussianFilter
% % % % % % % loglog(kx(1:Nkx).*2.*pi./Nx,E_fil(10,1:Nkx),'--k')
% % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % Spectral_GaussianFilter
% % % % % % % loglog(kx(1:Nkx).*2.*pi./Nx,EG(10,1:Nkx),'--r')
% % % % % % % xline(2.*pi/l);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % end
% % % % % % % if(0==1)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some matrix to zero
Flusso_Energia = zeros (Ny,Nx);
Flusso_Enstrofia = zeros (Ny,Nx);
% Derivatives set to zeros. 
Du_Dx  = zeros (Ny,Nx);
Dv_Dx  = zeros (Ny,Nx);
Du_Dy = zeros (Ny,Nx);
Dv_Dy = zeros (Ny,Nx);
DVort_Dx = zeros (Ny,Nx);
DVort_Dy = zeros (Ny,Nx);
% s count the interation on l.
s=0;
for idl=1:length(vl)
l=vl(idl);
s = s+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Differencies between Spatial & Spectral 
%                               &
%                       limit of our approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here I can choose between a spatial gaussian filter or a spectral filter
% that lead to results that are very similar but some differencies exist.
% Spectral filter appear to be much more efficient for cumputation rapidity
% as it goes much faster than the covolution product of spatial filters.
% It is also much more accurate at large scale when spatial filter do not
% manage to filter out the largest scales. Note also that spatial filter
% need to sistematically adjust the horizontal grid on which is computed
% the gaussian G(r), where r is the matrix size of the gaussian that the
% velocity field has to convoluate with. Spatial filter however is much
% more efficient to filter the smallest scales. When the velocity field is
% non periodic both technics consider strong velocity gradient at the edge 
% of the domain are small scales that have to be filtered out. This
% introduce very strong fluxes when the spatial derivative is computed in
% those outermost region. Consequently, in non periodic field, as it is the
% case for Torino data, we have to average the fluxes in the innermost par
% of the domain and cut out the edge. I see no other solution for that.
if(FSpatial==1) % set to 0 to have spatial filter and 1 to have spectral filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               Spatial Filter:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u_fil,v_fil,E_fil,Spakx,Spaky,SpaNkx,SpaNky] = func_Spatial_GaussianFilter(l,Nx,Ny,U,V,gaussian);
% if(0==1)
% loglog(kx(1:SpaNkx).*2.*pi./Nx,E_fil(10,1:SpaNkx),'--k')
% xline(2.*pi/l);
% end
[uu_fil,vv_fil,EE_fil,Spakx,Spaky,SpaNkx,SpaNky] = func_Spatial_GaussianFilter(l,Nx,Ny,U.*U,V.*V,gaussian); % U*U(y,x)
[vu_fil,uv_fil,ee_fil,Spakx,Spaky,SpaNkx,SpaNky] = func_Spatial_GaussianFilter(l,Nx,Ny,V.*U,U.*V,gaussian); % U*V(y,x)
% Vorticita
[Vort_fil,vort_fil,VorteeG,Spakx,Spaky,SpaNkx,SpaNky] = func_Spectral_GaussianFilter(l,Nx,Ny,vort,vort,gaussian); % U*V(y,x)
[Vortu_fil,Vortv_fil,VorteeG,Spakx,Spaky,SpaNkx,SpaNky] = func_Spectral_GaussianFilter(l,Nx,Ny,vort.*U,vort.*V,gaussian); % U*V(y,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              Spectral Filter:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
[u_fil,v_fil,EG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,U,V,gaussian); %U(y,x)
% if(0==1)
% loglog(kx(1:Nkx).*2.*pi./Nx,EG(10,1:Nkx),'--r')
% xline(2.*pi/l);
% end
[uu_fil,vv_fil,EEG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,U.*U,V.*V,gaussian); % U*U(y,x)
[vu_fil,uv_fil,eeG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,V.*U,U.*V,gaussian); % U*V(y,x)
% Vorticita
[Vort_fil,vort_fil,vorteeG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,vort,vort,gaussian); % U*V(y,x)
[Vortu_fil,Vortv_fil,vorteeG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,vort.*U,vort.*V,gaussian); % U*V(y,x)
end
% Imaginary part being zero
u_fil = real(u_fil);
v_fil = real(v_fil);
uu_fil = real(uu_fil);
vv_fil = real(vv_fil);
uv_fil = real(uv_fil);
vu_fil = real(vu_fil);
Vort_fil = real(Vort_fil);
Vortu_fil = real(Vortu_fil);
Vortv_fil = real(Vortv_fil);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                          DERIVATE D/Dx:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% ciclo per la derivata rispetto ad x
        for q = 1:Nx-1
            
            if q == 1
                
                Du_Dx (:,q) = (u_fil(:,q+1)-u_fil(:,q))/dx;
                Dv_Dx (:,q) = (v_fil(:,q+1)-v_fil(:,q))/dx;
                %
                DVort_Dx (:,q) = (Vort_fil(:,q+1)-Vort_fil(:,q))/dx;

               
            else
                
            Du_Dx (:,q) = (u_fil(:,q+1)-u_fil(:,q-1))/(2*dx);
            Dv_Dx (:,q) = (v_fil(:,q+1)-v_fil(:,q-1))/(2*dx);
            %
            DVort_Dx (:,q) = (Vort_fil(:,q+1)-Vort_fil(:,q-1))/(2*dx);

            
            end
            
        end
        
        for q = Nx
            
              Du_Dx (:,q) = (u_fil(:,q)-u_fil(:,q-1))/dx;
              Dv_Dx (:,q) = (v_fil(:,q)-v_fil(:,q-1))/dx;
              %
              DVort_Dx (:,q) = (Vort_fil(:,q)-Vort_fil(:,q-1))/dx;
            
        end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                          DERIVATE D/Dy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% ciclo per la derivata rispetto a y..dx=dy
        for p = 1: Ny-1
            
            if p == 1
                
                Du_Dy (p,:) = (u_fil (p+1,:)-u_fil (p,:))/dy;
                Dv_Dy (p,:) = (v_fil (p+1,:)-v_fil (p,:))/dy;
                %
                DVort_Dy (p,:) = (Vort_fil (p+1,:)-Vort_fil (p,:))/dy;
               
            else
                
            Du_Dy (p,:) = (u_fil (p+1,:)-u_fil (p-1,:))/(2*dy);
            Dv_Dy (p,:) = (v_fil (p+1,:)-v_fil (p-1,:))/(2*dy);
            %
            DVort_Dy (p,:) = (Vort_fil (p+1,:)-Vort_fil (p-1,:))/(2*dy);
            
            end
            
        end
        
        for p = Ny
            
              Du_Dy (p,:) = (u_fil (p,:)-u_fil (p-1,:))/dy;
              Dv_Dy (p,:) = (v_fil (p,:)-v_fil (p-1,:))/dy;
              %
              DVort_Dy (p,:) = (Vort_fil (p,:)-Vort_fil (p-1,:))/dy;
            
        end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                 FLUSSI:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%:.............................................................
        %%% Inicio ciclo per calcolare il flusso
        
        
        for i = 1: Ny
            for j = 1:Nx
                
                Flusso_Energia (i,j) = -(Du_Dx(i,j)*(uu_fil(i,j)-(u_fil(i,j)*u_fil(i,j)))+...
                                         Dv_Dy(i,j)*(vv_fil(i,j)-(v_fil(i,j)*v_fil(i,j)))+...
                                        (Du_Dy(i,j)+Dv_Dx(i,j))*(uv_fil(i,j)-(u_fil(i,j)*v_fil(i,j))));
                
                Flusso_Enstrofia (i,j) = -(DVort_Dx(i,j)*(Vortu_fil(i,j)-(u_fil(i,j)*Vort_fil(i,j)))+... 
                                           DVort_Dy(i,j)*(Vortv_fil(i,j)-(v_fil(i,j)*Vort_fil(i,j)))); 
                             
            end
        end
        
Flusso_l (:,t,s) = Flusso_Energia(:);
Flusso_Enstro_l (:,t,s) = Flusso_Enstrofia(:);
end
end
% Average in time.
Flusso_l_mt = squeeze(mean(Flusso_l,2));
Flusso_Enstro_l_mt = squeeze(mean(Flusso_Enstro_l,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idtc=60; % This is the truncation parameter.
% xmin=80;
% xmax=180;
% ymin=40;
% ymax=85;
% xmin=idxtc;
% xmax=Nx-idxtc;
Nx_tc=length([xmin:xmax]);
% ymin=idytc;
% ymax=Ny-idytc;
Ny_tc=length([ymin:ymax]);
% Energia
FLUSSO_Energia_tc = zeros (Ny_tc,Nx_tc,length(vl));
FLUSSO_Energia_l = zeros (length(vl),1);
FLUSSO_Energia_tl = zeros (length(vl),nFrames);
% Enstrofia
FLUSSO_Enstrofia_tc = zeros (Ny_tc,Nx_tc,length(vl));
FLUSSO_Enstrofia_l = zeros (length(vl),1);
FLUSSO_Enstrofia_tl = zeros (length(vl),nFrames);
for ids=1:length(vl)
    % Energia
    FLUSSO_Energia = reshape(Flusso_l_mt(:,ids),Ny,Nx);
    FLUSSO_Energia_tc(:,:,ids) = FLUSSO_Energia(ymin:ymax,xmin:xmax);
    FLUSSO_Energia_l(ids) = mean(mean(FLUSSO_Energia_tc(:,:,ids)));
    % Enstrofia
    FLUSSO_Enstrofia = reshape(Flusso_Enstro_l_mt(:,ids),Ny,Nx);
    FLUSSO_Enstrofia_tc(:,:,ids) = FLUSSO_Enstrofia(ymin:ymax,xmin:xmax);
    FLUSSO_Enstrofia_l(ids) = mean(mean(FLUSSO_Enstrofia_tc(:,:,ids)));
    for t=1:nFrames
        % Energia
        FLUSSO_Energia_t = reshape(Flusso_l(:,t,ids),Ny,Nx);
        FLUSSO_Energia_tl(ids,t) = mean(mean(FLUSSO_Energia_t(ymin:ymax,xmin:xmax)));
        % Enstrofia
        FLUSSO_Enstrofia_t = reshape(Flusso_Enstro_l(:,t,ids),Ny,Nx);
        FLUSSO_Enstrofia_tl(ids,t) = mean(mean(FLUSSO_Enstrofia_t(ymin:ymax,xmin:xmax)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    PLOTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we plot the total field on which the filtering is performed and the
% truncated field on which we averaged the fluxes, i.e. getting rid off the
% edge of the domain.
if(0==0)
figure
contourf(X(ymin:ymax,xmin:xmax),Y(ymin:ymax,xmin:xmax),FLUSSO_Energia_tc(:,:,120), 250, 'LineStyle','none');
title('Energy flux')
colorbar
figure
contourf(X,Y,U, 250, 'LineStyle','none');
hold on
scatter(X(ymin,xmin),Y(ymin,xmin), [], 'Xr')
scatter(X(ymin,xmax),Y(ymin,xmax), [], 'Xr')
scatter(X(ymax,xmin),Y(ymax,xmin), [], 'Xr')
scatter(X(ymax,xmax),Y(ymax,xmax), [], 'Xr')
colorbar
title('Full field - crosses are truncation')
figure
contourf(X,Y,u_fil, 250, 'LineStyle','none');
colorbar
title('U filterd')
figure
contourf(X(ymin:ymax,xmin:xmax),Y(ymin:ymax,xmin:xmax),U(ymin:ymax,xmin:xmax), 250, 'LineStyle','none');
colorbar
title('Truncated field')
end
%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   SAVE
%  ###########################################################################################################################
%  ###########################################################################################################################
% FLUSSO_Energia_l(length(vl),1)
% FLUSSO_Energia_tl(length(vl),t)
Fluxes = [FLUSSO_Energia_l,FLUSSO_Enstrofia_l];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time averaged
fileout = [roots,Name,'/Fluxes_mt'] ;
filename = sprintf('%s.mtx',fileout);
fid = fopen(filename,'wb');
fwrite(fid,size(Fluxes,1),'ulong');
fwrite(fid,size(Fluxes,2),'ulong');
fwrite(fid,Fluxes(:),'float');
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time dependant
% % % 
fileout1 = [roots,Name,'/Energy_flux_Fr_',num2str(nTime)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(FLUSSO_Energia_tl,1),'ulong');
fwrite(fid,size(FLUSSO_Energia_tl,2),'ulong');
fwrite(fid,FLUSSO_Energia_tl(:),'float');
fclose(fid);
% % % 
fileout1 = [roots,Name,'/Enstrophy_flux_Fr_',num2str(nTime)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(FLUSSO_Enstrofia_tl,1),'ulong');
fwrite(fid,size(FLUSSO_Enstrofia_tl,2),'ulong');
fwrite(fid,FLUSSO_Enstrofia_tl(:),'float');
fclose(fid);
% % % 
save([roots,Name,'/EFluxes_Spectral_infos.mat'],'Nkx','Nky','nFrames','x','y','kx','ky','Nx','Ny','vl','dx','dy')