%% script per plottare il flusso di energia e quello di enstrofia
%% considerando la media spaziale solo nel settore circolare considerato, 
%% singolarmente e sovrapposti
nodi =128;
   
dth = 360/Nraggi;

%calcolo il numero del raggio corrispondente ad un dato angolo 
R1 = round(Ntmin/dth);
R2 = round(Ntmax/dth);

name_mtx = [roots,Name,'/Energy_flux_dl_' num2str(dl)];
[F_Energy1,b] = loadmtx (name_mtx);

name_mtx2 = [roots,Name,'/Enstrophy_flux_dl_' num2str(dl)];
[F_Enstrophy,b] = loadmtx (name_mtx2);

%creo griglia cartesiana...verr� una griglia di 128*128 punti con gli
%estremi giusti
[X,Y]=Create_Grid_cart_Vasca('C',[nodi,nodi,-r_max_cm,r_max_cm,-r_max_cm,r_max_cm]);


L_max = 2*r_max_cm;
N_l = round(4/dl + (L_max-4));

asse_l = [(dl:dl:4) 5:L_max] ;
asse_k = (2*pi)./asse_l ;

totaltime = b(1)/N_l/(nodi*nodi);

%----------------------------------------------------------------------
%reshape flusso di energia
%----------------------------------------------------------------------
F_Energy = reshape(F_Energy1,nodi*nodi,round(totaltime),N_l);

%F_Energy2 = reshape(F_Energy2,nodi*nodi,round(totaltime),N_l);

%F_Energy =[F_Energy1 F_Energy2];

%media temporale
Fm=mean(F_Energy,2);

%ciclo di interpolazione di Fm per fare media nel settore circolare
Fm_n_pol_m = zeros (1,N_l);

for n=1:N_l
    
    Fm_n = reshape (Fm (:,1,n),nodi,nodi);
    
    % interpolazione su griglia polare del flusso per un dato n
    Fm_n_pol = interp2 ( X, Y,Fm_n,Grid_Xp_cm,Grid_Yp_cm,'spline');
    
    %media spaziale nel settore circolare considerato
    Fm_n_pol_m (n) = mean (mean(Fm_n_pol (R1:R2,:)));
  
end   
    
%%---------------------------------------------------------------------
figure  %Energy & asse k

  scrsz = get(0,'ScreenSize');

   set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')

   axes('FontSize',18,'Linewidth',2,'FontName','times',...
       'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
 
   hold on
 
   subplot ('Position',[0.13 0.17 0.82 0.75]) 
  
    
   hold on
plot(asse_k,Fm_n_pol_m,'k','linewidth',1)
%plot(kz, Ez_k,'r','linewidth',1)
   
 %axis([10^(-1) 10^2 10^(-5) 10^0 ])

ylabel ('$\Pi$(cm$^2$s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')

 axis ([10^-1 10^2 -0.002 0.002])

 
 %% ---------------------------------------------------------------------
 figure   % Energy & asse l

  scrsz = get(0,'ScreenSize');

   set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')

   axes('FontSize',18,'Linewidth',2,'FontName','times',...
       'TickLength',[0.01; 0.03],'Xscale','Lin',...
     'Position',[0.13 0.17 0.82 0.75])   
 
   hold on
 
   subplot ('Position',[0.13 0.17 0.82 0.75]) 
  
    
   hold on
plot(asse_l,Fm_n_pol_m,'k','linewidth',1)
%plot(kz, Ez_k,'r','linewidth',1)
   
 %axis([10^(-1) 10^2 10^(-5) 10^0 ])

ylabel ('$\Pi$(cm$^2$s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$l$(cm)','FontSize',18,'FontName','times','Interpreter','Latex')

 axis ([0.4 60 -0.002 0.002])
 
%-----------------------------------------------------------------------
%reshape flusso di enstrofia
%-----------------------------------------------------------------------
F_Enstrophy = reshape(F_Enstrophy,nodi*nodi,totaltime,N_l);

%media spaziale
%Fm_2=mean(F_Enstrophy,1);

%media temporale
Fm2=mean(F_Enstrophy,2);

%flusso di enstrofia in funzione delle scale di cut-off
%Fm_2=Fm_2(:);

%ciclo di interpolazione di Fm per fare media nel settore circolare
Hm_n_pol_m = zeros (1,N_l);

for n=1:N_l
    
    Fm2_n = reshape (Fm2 (:,1,n),nodi,nodi);
    
    % interpolazione su griglia polare del flusso per un dato n
    Fm2_n_pol = interp2 ( X, Y,Fm2_n,Grid_Xp_cm,Grid_Yp_cm,'spline');
    
    %media spaziale nel settore circolare considerato
    Hm_n_pol_m (n) = mean (mean(Fm2_n_pol (R1:R2,:)));
  
end  

%%-----------------------------------------------------------------------
figure  %% Enstrophy & asse k

  scrsz = get(0,'ScreenSize');

   set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')

   axes('FontSize',18,'Linewidth',2,'FontName','times',...
       'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
 
   hold on
 
   subplot ('Position',[0.13 0.17 0.82 0.75]) 
  
    
   hold on
plot(asse_k,Hm_n_pol_m,'k','linewidth',1)
%plot(kz, Ez_k,'r','linewidth',1)
   

 axis([4*10^(-1) 10^2 -0.0002 0.0002 ])

ylabel ('$Z$(s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')

%%-----------------------------------------------------------------------
figure  %% Enstrophy & asse l

  scrsz = get(0,'ScreenSize');

   set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')

   axes('FontSize',18,'Linewidth',2,'FontName','times',...
       'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
 
   hold on
 
   subplot ('Position',[0.13 0.17 0.82 0.75]) 
  
    
   hold on
plot(asse_l,Hm_n_pol_m,'k','linewidth',1)
%plot(kz, Ez_k,'r','linewidth',1)
   

 axis([4*10^(-1) 10^2 -0.0002 0.0002 ])

ylabel ('$Z$(s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$l$(cm)','FontSize',18,'FontName','times','Interpreter','Latex')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    FLUSSI SOPRAPOPOSTI:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%--------------------flussi sovrapposti in l --------------------------------

figure

  scrsz = get(0,'ScreenSize');

   set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')

   axes('FontSize',18,'Linewidth',2,'FontName','times',...
       'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
 
   hold on
 
   subplot ('Position',[0.13 0.17 0.82 0.75]) 
  
    
   hold on

%plot(asse_k,Fm_n_pol_m,'-*k','linewidth',1)
%plot(asse_k,Fm2_n_pol_m,'-sk','linewidth',1)

plot(asse_l,Fm_n_pol_m,'-*r','linewidth',1.5)

plot(asse_l,Hm_n_pol_m,'-sb','linewidth',1.5)
%plot(kz, Ez_k,'r','linewidth',1)
   

 axis([4*10^(-1) 10^2 -0.0015 0.0015 ])

ylabel ('$\Pi$(cm$^2$s$^{-3}$), $Z$(s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
%xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$l$(cm)','FontSize',18,'FontName','times','Interpreter','Latex')


%%%--------------------flussi sovrapposti in k --------------------------------

figure

  scrsz = get(0,'ScreenSize');

   set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')

   axes('FontSize',18,'Linewidth',2,'FontName','times',...
       'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
 
   hold on
 
   subplot ('Position',[0.13 0.17 0.82 0.75]) 
  
    
   hold on

%plot(asse_k,Fm_n_pol_m,'-*k','linewidth',1)
%plot(asse_k,Fm2_n_pol_m,'-sk','linewidth',1)

plot(asse_k,Fm_n_pol_m,'-*r','linewidth',1.5)

plot(asse_k,Hm_n_pol_m,'-sb','linewidth',1.5)
%plot(kz, Ez_k,'r','linewidth',1)
   

 axis([asse_k(end) asse_k(1) -0.0015 0.0015 ])

ylabel ('$\Pi$(cm$^2$s$^{-3}$), $Z$(s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
%xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')