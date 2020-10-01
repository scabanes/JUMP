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
%---------------------------------------------
% Here I used Roland's interpolation through the data merged.civ2. The
% following gapt allow to load data in time.
% % % % gapt = 30;
criterionNaN = 20; % UNSUED -- window on which the mean to replace Nan is calculated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % ---------------------------------------------------- for pseudo-cartesian Torino
% % % % % % % % % % % % % % % % % % % % % % % % if(1==0)
% % % % % % % % % % % % % % % % % % % % % % % % % Here we load velocity fields:
% % % % % % % % % % % % % % % % % % % % % % % % [Vtheta_tot,b]=loadmtx([roots,Name,NameVt]);
% % % % % % % % % % % % % % % % % % % % % % % % [Vr_tot,b]=loadmtx([roots,Name,NameVr]);
% % % % % % % % % % % % % % % % % % % % % % % % [curl_tot,b]=loadmtx([roots,'EXPT03','/VortRelat_time_540_180_1105']);
% % % % % % % % % % % % % % % % % % % % % % % % % Here we map the loaded fields:
% % % % % % % % % % % % % % % % % % % % % % % % [returnOK] = Maps(Vtheta_tot,Vr_tot,Grid_Xp_cm,Grid_Yp_cm,GridR_2C,GridT_2C,Nti,Nri,Nrmin,Nrmax,Ntmin,Ntmax,theta,1)
% % % % % % % % % % % % % % % % % % % % % % % % % We truncated the domain size if necessary.
% % % % % % % % % % % % % % % % % % % % % % % % Grid_Xp_cm=Grid_Xp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
% % % % % % % % % % % % % % % % % % % % % % % % Grid_Yp_cm=Grid_Yp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
% % % % % % % % % % % % % % % % % % % % % % % % % Locally pseudo-cartesian grid see Read. 2015
% % % % % % % % % % % % % % % % % % % % % % % % GridR_2C=GridR_2C(Ntmin:Ntmax,Nrmin:Nrmax);
% % % % % % % % % % % % % % % % % % % % % % % % GridT_2C=GridT_2C(Ntmin:Ntmax,Nrmin:Nrmax);
% % % % % % % % % % % % % % % % % % % % % % % % % -----------------------------------------------------Truncation
% % % % % % % % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % % % % % % % % ---------------------------------------------------- for pseudo-cartesian Torino
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
% Flusso_l = zeros (length([idymin:idymax])*length([idxmin:idxmax]),nFrames,length(vl));
Flusso_l = zeros (length([idymin:idymax])*length([idxmin:idxmax]),[],length(vl));
Flusso_Enstro_l = zeros (length([idymin:idymax])*length([idxmin:idxmax]),[],length(vl));
% ##################################################################################################################################################
%                                                                                                                                   INTERPOLAZIONE :
% ##################################################################################################################################################
it=0;
for t=Itime:Tmax %---------------------------------------------------------------------------Init: LOOP on t
    
% % % % % % % % if(0==0)%------------------------------------------------------- Cas Torino
% file to load 
% file = [roots,Name,'/merged.civ2/fig_',num2str(numt),'.nc'];
% if isfile(file)
%      numt = numt+gapt;% File exists.
% % compter les images
% % Velocity fields of cartesian velocities
% UU = double(ncread(file,'U'));
% VV = double(ncread(file,'V'));
% curl = double(ncread(file,'curl'));
% % regular cartesian coordinates
% xx = double(ncread(file,'coord_x'));
% yy = double(ncread(file,'coord_y'));
% 
% x = xx(idxmin:idxmax);
% y = yy(idymin:idymax);
% % % % % % % % % % % % % % % %----------------------------------------------for pseudo-cartesian Torino
% % % % % % % % % % % % % % % % % % % % % % % % % % UU = reshape(Vtheta_tot(:,t),Nti,Nri);
% % % % % % % % % % % % % % % % % % % % % % % % % % VV = reshape(Vr_tot(:,t),Nti,Nri);
% % % % % % % % % % % % % % % % % % % % % % % % % % curl = reshape(curl_tot(:,t),Nti,Nri);
% % % % % % % % % % % % % % % % % % % % % % % % % % %----------------- Truncation of the velocity fields
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % U=UU(Ntmin:Ntmax,Nrmin:Nrmax);%.*Tocm;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % V=VV(Ntmin:Ntmax,Nrmin:Nrmax);%.*Tocm;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % vort=curlL(Ntmin:Ntmax,Nrmin:Nrmax);%.*Tocm;
% % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % theta_tc = theta(1:length(idymin:idymax)); % Attenzione theta deve partire da zero..
% % % % % % % % % % % % % % % % % % % % % % % % % % r_tc = r(idxmin:idxmax); % Attenzione theta deve partire da zero..
% % % % % % % % % % % % % % % % % % % % % % % % % % % Here we choose a pseudo cartesian grid 
% % % % % % % % % % % % % % % % % % % % % % % % % % xx = Ro.*theta; 
% % % % % % % % % % % % % % % % % % % % % % % % % % yy = r;
% % % % % % % % % % % % % % % % % % % % % % % % % % x = Ro.*theta_tc;
% % % % % % % % % % % % % % % % % % % % % % % % % % y = r_tc;
% % % % % % % % % % % % % % % %--------------------------------------------for pseudo-cartesian Torino
% 
% Nx = length(x);
% Ny = length(y);
% dx=mean(x(2:end)-x(1:end-1));
% dy=mean(y(2:end)-y(1:end-1));
% Lx = x(end)-x(1); % en m
% Ly = y(end)-y(1); % en m
% % Recreat the mesh grid
% [XX,YY] = meshgrid(xx,yy);

% % % % % % % % % % elseif(1==0)%------------------------------------------------------- Cas Obs Jupiter
% % % % % % % % % % % file to load    
% % % % % % % % % % file = [roots,Name,'StatisticalData.nc'];
% % % % % % % % % % % Velocity fields of cartesian velocities
% % % % % % % % % % UUL = double(ncread(file,'wm')); %wm(lat,long,z,time)
% % % % % % % % % % VVL = double(ncread(file,'vm')); %vm(lat,long,z,time)
% % % % % % % % % % curlL = double(ncread(file,'vort')); %vort(lat,long,z,time)
% % % % % % % % % % % single time
% % % % % % % % % % UU = UUL(:,:,1,t);
% % % % % % % % % % VV = VVL(:,:,1,t);
% % % % % % % % % % curl = curlL(:,:,1,t);
% % % % % % % % % % %
% % % % % % % % % % SU=size(UU);
% % % % % % % % % % % regular cartesian coordinates
% % % % % % % % % % xx=[1:(2.*pi.*R_Jup)./SU(2):2.*pi.*R_Jup];
% % % % % % % % % % yy=[1:(pi.*R_Jup)./SU(1):pi.*R_Jup];
% % % % % % % % % % %
% % % % % % % % % % x = xx(idxmin:idxmax);
% % % % % % % % % % y = yy(idymin:idymax);
% % % % % % % % % % Nx = length(x);
% % % % % % % % % % Ny = length(y);
% % % % % % % % % % dx=mean(x(2:end)-x(1:end-1));
% % % % % % % % % % dy=mean(y(2:end)-y(1:end-1));
% % % % % % % % % % Lx = x(end)-x(1); % en m
% % % % % % % % % % Ly = y(end)-y(1); % en m
% % % % % % % % % % % Recreat the mesh grid
% % % % % % % % % % [XX,YY] = meshgrid(xx,yy);
% % % % % % % % % % 
% % % % % % % % % % 
% % % % % % % % % % else%------------------------------------------------------- Cas Vasca
[UUL,b]=loadmtx([roots,Name,NameU]); % It has to be checked, but this is Uy
[VVL,b]=loadmtx([roots,Name,NameV]); % It has to be checked, but this is Uy
[curlL,b]=loadmtx([roots,Name,NameVort]);
% We work on a squared cartesian grid. nodi is the length of it. 
nodi = b(1)^0.5;
[XX,YY]=Create_Grid_cart_Vasca('C',[nodi,nodi,-r_max_cm,r_max_cm,-r_max_cm,r_max_cm]);
%
UU = reshape(UUL(:,t),nodi,nodi);
VV = reshape(VVL(:,t),nodi,nodi);
curl = reshape(curlL(:,t),nodi,nodi);
% % % % % % % % % % end

% Cut to a rectangular table
U = UU(idymin:idymax,idxmin:idxmax);
V = VV(idymin:idymax,idxmin:idxmax);
vort = curl(idymin:idymax,idxmin:idxmax);
X = XX(idymin:idymax,idxmin:idxmax);
Y = YY(idymin:idymax,idxmin:idxmax);
clear VV curl
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                    NAN VALUES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
it = it + 1;
% count Nan
NbNan_U(it)=length(find(isnan(U)==1));
% NbNan_V(t)=length(find(isnan(V)==1));
% NbNan_Vort(t)=length(find(isnan(vort)==1));
% Fill in the Nans in the table with 0 values
U(isnan(U))=0.;
V(isnan(V))=0.;
vort(isnan(vort))=0.;
% This is related to the saving of the fluxes at the end of the loops
% If toot much Nan values
% if(NbNan_U(t)>NbMaxNan)  %---------------------------------------------------------------------------Init: Condition on Nan
% it = it;
% disp(['Save frame ',num2str(it),' on frames ',num2str(t)])
% disp('Too much Nan values')
% else
disp(['Save frame ',num2str(it),' on frames ',num2str(t)])
% if(1==0)
%  figure
%  contourf(X,Y,U, 100, 'LineStyle','none');
%  colorbar
%  size(U)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        Field spectral decomposition:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        FOURIER DECOMPOSITION for U:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
loglog(kx(1:Nkx).*2.*pi./Nx,diag(E(1:Nky,1:Nkx)))
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
%                                                                                                            LOOP on length l:
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
for idl=1:length(vl) %---------------------------------------------------------------------------Init: LOOP on l
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
% % loglog(kx(1:SpaNkx).*2.*pi./Nx,E_fil(10,1:SpaNkx),'--k')
% % xline(2.*pi/l);
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
        if(l==1) % No filter for l=grid=1
            u_fil = U;
            v_fil = V;
            uu_fil = U.*U;
            vv_fil = V.*V;
            vu_fil =  V.*U;
            uv_fil = U.*V ;
            Vort_fil = vort;
            Vortu_fil = vort.*U;
            Vortv_fil = vort.*V;
        else
[u_fil,v_fil,EG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,U,V,gaussian); %U(y,x)
% if(0==1)
% loglog(kx(1:Nkx).*2.*pi./Nx,diag(EG(1:Nky,1:Nkx)),'--r')
% xline(2.*pi/l);
% end
[uu_fil,vv_fil,EEG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,U.*U,V.*V,gaussian); % U*U(y,x)
[vu_fil,uv_fil,eeG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,V.*U,U.*V,gaussian); % U*V(y,x)
% Vorticita
[Vort_fil,vort_fil,vorteeG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,vort,vort,gaussian); % U*V(y,x)
[Vortu_fil,Vortv_fil,vorteeG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,vort.*U,vort.*V,gaussian); % U*V(y,x)
        end
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

  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Save frames (without too much Nan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Flusso_l (:,it,s) = Flusso_Energia(:);
Flusso_Enstro_l (:,it,s) = Flusso_Enstrofia(:);


end%------------------------------------------------------------------------------------- End: LOOP on l
% end%------------------------------------------------------------------------------------- End: condition if too much Nan
clear Du_Dx Dv_Dx Du_Dy Dv_Dy DVort_Dx DVort_Dy 
% else% File does not exist.
%     numt = numt + gapt;
% end%------------------------------------------------------------------------------------- End: file exists
end%--------------------------------------------------------------------------------------End: LOOP on t
nbFinalt = it;
clear t
disp(['Number of frames:',num2str(nbFinalt)])
% Average in time.
Flusso_l_mt = squeeze(mean(Flusso_l,2));
Flusso_Enstro_l_mt = squeeze(mean(Flusso_Enstro_l,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               AVERAGE II:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==0)%----------------------------------- average 2 Options
% un critere supllementaire pour suprimer la partie interne de la zone a
% moyenner.
Ax=[xmin:xmax];
Ay=[ymin:ymax];
%
if(radius==0)
circlePixels = ones(Ny,Nx);
else
[columnsInImage rowsInImage] = meshgrid(1:Ny, 1:Nx);
% Next create the circle in the image.
% radius = 10;
centerX = Ny/2;
centerY = Nx/2;
circlePixels = double( (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2);
circlePixels(find(circlePixels==1))=nan;
circlePixels(find(circlePixels==0))=1;
end
%
Nx_tc=length(Ax);
Ny_tc=length(Ay);
%--------------------------------------------------------------------------
% Energia
FLUSSO_Energia_tc = zeros (Ny_tc,Nx_tc,length(vl));
FLUSSO_Energia_l = zeros (length(vl),1);
FLUSSO_Energia_tl = zeros (length(vl),[]);
% Enstrofia
FLUSSO_Enstrofia_tc = zeros (Ny_tc,Nx_tc,length(vl));
FLUSSO_Enstrofia_l = zeros (length(vl),1);
FLUSSO_Enstrofia_tl = zeros (length(vl),[]);
%--------------------------------------------------------------------------
for ids=1:length(vl)
    % Energia
    FLUSSO_Energia = reshape(Flusso_l_mt(:,ids),Ny,Nx);
    FLUSSO_Energia = FLUSSO_Energia.*circlePixels;
    FLUSSO_Energia_tc(:,:,ids) = FLUSSO_Energia(Ay,Ax);
    FLUSSO_Energia_NN = FLUSSO_Energia(Ay,Ax);
    NoNan = ~isnan(FLUSSO_Energia_NN);
    FLUSSO_Energia_l(ids) = mean(FLUSSO_Energia_NN(NoNan));
%     FLUSSO_Energia_l(ids) = mean(mean(FLUSSO_Energia_tc(:,:,ids)));
    % Enstrofia
    FLUSSO_Enstrofia = reshape(Flusso_Enstro_l_mt(:,ids),Ny,Nx);
    FLUSSO_Enstrofia = FLUSSO_Enstrofia.*circlePixels;
    FLUSSO_Enstrofia_tc(:,:,ids) = FLUSSO_Enstrofia(Ay,Ax);
    FLUSSO_Enstrofia_NN = FLUSSO_Enstrofia(Ay,Ax);
    NoNan = ~isnan(FLUSSO_Enstrofia_NN);
    FLUSSO_Enstrofia_l(ids) = mean(FLUSSO_Enstrofia_NN(NoNan));
%     FLUSSO_Enstrofia_l(ids) = mean(mean(FLUSSO_Enstrofia_tc(:,:,ids)));
    for t=1:nbFinalt
        % Energia
        FLUSSO_Energia_t = reshape(Flusso_l(:,t,ids),Ny,Nx);
        FLUSSO_Energia_t = FLUSSO_Energia_t.*circlePixels;
        FLUSSO_Energia_t = FLUSSO_Energia_t(Ay,Ax);
        NoNan = ~isnan(FLUSSO_Energia_t);
        FLUSSO_Energia_tl(ids,t) = mean(FLUSSO_Energia_t(NoNan));
%         FLUSSO_Energia_tl(ids,t) = mean(mean(FLUSSO_Energia_t(Ay,Ax)));
        % Enstrofia
        FLUSSO_Enstrofia_t = reshape(Flusso_Enstro_l(:,t,ids),Ny,Nx);
        FLUSSO_Enstrofia_t = FLUSSO_Enstrofia_t.*circlePixels;
        FLUSSO_Enstrofia_t = FLUSSO_Enstrofia_t(Ay,Ax);
        NoNan = ~isnan(FLUSSO_Enstrofia_t);
        FLUSSO_Energia_tl(ids,t) = mean(FLUSSO_Enstrofia_t(NoNan));
%         FLUSSO_Enstrofia_tl(ids,t) = mean(mean(FLUSSO_Enstrofia_t(Ay,Ax)));
    end
end
end%----------------------------------- End average 2 Options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    PLOTS: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we plot the total field on which the filtering is performed and the
% truncated field on which we averaged the fluxes, i.e. getting rid off the
% edge of the domain.
if(0==0)
figure
contourf(XX,YY,UU, 250, 'LineStyle','none');
hold on
scatter(XX(idymin,idxmin),YY(idymin,idxmin), [], 'Xr')
scatter(XX(idymin,idxmax),YY(idymin,idxmax), [], 'Xr')
scatter(XX(idymax,idxmin),YY(idymax,idxmin), [], 'Xr')
scatter(XX(idymax,idxmax),YY(idymax,idxmax), [], 'Xr')
colorbar
title('Full field - crosses are truncation')
%
figure
contourf(X(Ay,Ax),Y(Ay,Ax),FLUSSO_Energia_tc(:,:,2), 250, 'LineStyle','none');
% imagesc(FLUSSO_Energia_tc(:,:,2));
title('Energy flux - small scale in the averaged surface')
colorbar
caxis([-0.01 0.01]);
%
figure
contourf(X(Ay,Ax),Y(Ay,Ax),FLUSSO_Energia_tc(:,:,7), 250, 'LineStyle','none');
title('Energy flux - medium scale in the averaged surface')
colorbar
caxis([-0.01 0.01]);
%
figure
contourf(X(Ay,Ax),Y(Ay,Ax),FLUSSO_Energia_tc(:,:,17), 250, 'LineStyle','none');
title('Energy flux - large scale  in the averaged surface')
colorbar
caxis([-0.01 0.01]);
%
figure
contourf(X,Y,U, 250, 'LineStyle','none');
hold on
scatter(X(ymin,xmin),Y(ymin,xmin), [], 'Xr')
scatter(X(ymin,xmax),Y(ymin,xmax), [], 'Xr')
scatter(X(ymax,xmin),Y(ymax,xmin), [], 'Xr')
scatter(X(ymax,xmax),Y(ymax,xmax), [], 'Xr')
colorbar
title('Full truncated field - crosses are averaging truncation')
% figure
% contourf(X,Y,u_fil, 250, 'LineStyle','none');
% colorbar
% title('U filterd')
figure
contourf(X(Ay,Ax),Y(Ay,Ax),U(Ay,Ax), 250, 'LineStyle','none');
colorbar
title('Averaging truncated field')
% nombre de Nan
figure; hold on
plot(NbNan_U,'linewidth',2)
% plot(NbNan_V,'linewidth',2)
% plot(NbNan_Vort,'linewidth',2)
legend('Nan on U')

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
fileout = [roots,Name,'/Fluxes_mt_Fr_',num2str(nTime),'_Itime_',num2str(Itime),'_Idymin_',num2str(idymin)] ;
filename = sprintf('%s.mtx',fileout);
fid = fopen(filename,'wb');
fwrite(fid,size(Fluxes,1),'ulong');
fwrite(fid,size(Fluxes,2),'ulong');
fwrite(fid,Fluxes(:),'float');
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time dependant
% % % 
fileout1 = [roots,Name,'/Energy_flux_Fr_',num2str(nTime),'_Itime_',num2str(Itime),'_Idymin_',num2str(idymin)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(FLUSSO_Energia_tl,1),'ulong');
fwrite(fid,size(FLUSSO_Energia_tl,2),'ulong');
fwrite(fid,FLUSSO_Energia_tl(:),'float');
fclose(fid);
% % % 
fileout1 = [roots,Name,'/Enstrophy_flux_Fr_',num2str(nTime),'_Itime_',num2str(Itime),'_Idymin_',num2str(idymin)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(FLUSSO_Enstrofia_tl,1),'ulong');
fwrite(fid,size(FLUSSO_Enstrofia_tl,2),'ulong');
fwrite(fid,FLUSSO_Enstrofia_tl(:),'float');
fclose(fid);
% % % 
save([roots,Name,'/EFluxes_Spectral_infos.mat'],'Nkx','Nky','nFrames','Itime','kx','ky','Nx','Ny','vl','dx','dy')