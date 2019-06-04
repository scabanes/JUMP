%clear all

% Create_Grid_pol_OK
% nodi =128;
% 
% %creo griglia cartesiana...verrà una griglia di 128*128 punti con gli
% %estremi giusti
% [X,Y]=Create_Grid_cart_pol_OK('C',[nodi,nodi,-r_max_cm,r_max_cm,-r_max_cm,r_max_cm]);


 salva_movie=0;      % si=1 no=0
 if(salva_movie==1)
%      mov = avifile('VortRelat_movie_ESP9.avi');
%       
%      mov.compression= 'None';
%      mov.Quality = 100;
%      mov.fps=20;

      mov=VideoWriter('VortPot_movie_ESP9.avi');
      mov.Quality = 100;
      open(mov)
 end


% n_in = input('Number of the mtr(1, 2...)?: ','s');
% n =sscanf(n_in,'%f');

name_mtx = ['VortPot_time_04'];
[vpot,b] = loadmtx (name_mtx); 
[U,b] = loadmtx ('U5_pol_540_160_2041'); 
[V,b] = loadmtx ('V5_pol_540_160_2041'); 
 
%se ci vogliamo sovrapporre i files polari?
% name_mtx1 = ['Vz'];
% [Vz,b] = loadmtx (name_mtx1);
%[Vz,b] = loadmtx (name_mtx1); 
% name_mtx2 = ['Vr'];
% [Vr,b] = loadmtx (name_mtx2);
% name_mtx2 = ['Vr'];

% Load all files
 %o raw?
% U_out = double(ncread([path_out,'fig_',num2str(t),'-',num2str(t+1),'.nc'],'Civ2_U'));
% V_out = double(ncread([path_out,'fig_',num2str(t),'-',num2str(t+1),'.nc'],'Civ2_V'));
% X_out = double(ncread([path_out,'fig_',num2str(t),'-',num2str(t+1),'.nc'],'Civ2_X'));
% Y_out = double(ncread([path_out,'fig_',num2str(t),'-',num2str(t+1),'.nc'],'Civ2_Y'));
% U_in = double(ncread([path_in,'fig_',num2str(t),'-',num2str(t+1),'.nc'],'Civ2_U'));
% V_in = double(ncread([path_in,'fig_',num2str(t),'-',num2str(t+1),'.nc'],'Civ2_V'));
% X_in = double(ncread([path_in,'fig_',num2str(t),'-',num2str(t+1),'.nc'],'Civ2_X'));
% Y_in = double(ncread([path_in,'fig_',num2str(t),'-',num2str(t+1),'.nc'],'Civ2_Y'));
% % concatenate in one vector inner and outer camera
% U = [U_out;U_in];
% V = [V_out;V_in];
% X = [X_out;X_in];
% Y = [Y_out;Y_in];
% name_mtx1 = ['Vz'];

  
  step = 1;
  t=0;
  
      [n,T]=size(vpot);
      totaltime = T;
      
      for i = 1:step:totaltime
          
          t=t+1
          vpot_i=reshape(vpot(:,i),Nraggi,Ncerchi);
          %max=max(abs((vrel_i(:))));
        
%           Vz_i=reshape(U(:,i),nodi,nodi);
%           Vr_i=reshape(V(:,i),nodi,nodi);
   
       figure(1), contourf(Grid_Xp_m,Grid_Yp_m,vpot_i,250,'LineStyle','none')
% % %        set(gca,'XTick',[],'YTick',[],'LineWidth',1,...
% % %       'DataAspectRatio',[1 1 1])%,'Xticklabel',[])
        axis([-1 1.5 0.5 2.5]);
      colormap([fliplr(hot(40));  ])
      colorbar ('Fontsize',14,'FontName','times')%,'YTickLabel',([-1 0 1]))
      caxis([-0.5 0.8])
      hold on
      quiver(Grid_Xp_m,Grid_Yp_m,reshape(U(:,i),540,160),reshape(V(:,i),540,160),'k');
   % pause
   %hold on
    note = annotation('textbox',[0.77 0.9 0.1 0.1]);
   %set(note,'String','${\psi$(cm$^$2s$^{-1})}$ ','Interpreter','Latex','FontSize',16,'LineStyle','none') 
   set(note,'String','${vort_pot}$ ','Interpreter','Latex','FontSize',16,'LineStyle','none') 
   pause
%    hold on
%    
%    %streamslice(X,Y,U_i,V_i,4)
%    quiver (X(1:2:end,1:2:end),Y(1:2:end,1:2:end),U_i(1:2:end,1:2:end),V_i(1:2:end,1:2:end),2,'k')
   
   
   % title(sprintf('Tempo= %f min', tempo/60))
    
    if(salva_movie==1)
    
            F = getframe(gcf);
            %filename=sprintf('%s\\%s%06d.tif','c:','out',k);
            %filename=[fig num2str(i,'%.5i')];
            %saveas(gcf, filename, 'jpeg')
            %mov = addframe(mov,F);
            writeVideo(mov,F);

            
    end
    
    %pause
    figure(1),clf
    clear max
      end
 
  
if(salva_movie==1)    
    close(mov);
end

 fclose all
 

        
        
        
        
        
        