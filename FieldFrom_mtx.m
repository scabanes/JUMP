%########################################################################################################
%                                                                                           LOAD it -> 0:
%########################################################################################################
if(it==0)
[ux,b]=loadmtx([roots,Name,NameU]); % It has to be checked, but this is Ux
[uy,b]=loadmtx([roots,Name,NameV]); % It has to be checked, but this is Uy
[curlL,b]=loadmtx([roots,Name,NameVort]);
% We work on a squared cartesian grid. nodi is the length of it. 
nodi = b(1)^0.5;
[XX,YY]=Create_Grid_cart_Vasca('C',[nodi,nodi,-r_max_cm,r_max_cm,-r_max_cm,r_max_cm]);
%
Ux = reshape(ux(:,1),nodi,nodi);
Uy = reshape(uy(:,1),nodi,nodi);%-
curl = reshape(curlL(:,1),nodi,nodi);
[returnOK] = Maps_Cartesian(Ux,Uy,XX,YY,idxmin,idxmax,idymin,idymax)
end
%########################################################################################################
%                                                                                           LOAD it -> 1:
%########################################################################################################
if(it>0)
[ux,b]=loadmtx([roots,Name,NameU]); % It has to be checked, but this is Ux
[uy,b]=loadmtx([roots,Name,NameV]); % It has to be checked, but this is Uy
[curlL,b]=loadmtx([roots,Name,NameVort]);
% We work on a squared cartesian grid. nodi is the length of it. 
nodi = b(1)^0.5;
%
Ux = reshape(ux(:,it),nodi,nodi);
Uy = reshape(uy(:,it),nodi,nodi);%-
curl = reshape(curlL(:,it),nodi,nodi);
% % % % % % % % % % end
end