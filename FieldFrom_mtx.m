%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mostly when derivative are involved, e.g energy and enstrophy fluxes.
% It is of prime interest to well define the sign of the velocity field on
% axis x and y as well as the direction of the derivatives along this axis
% by setting 'DevType_x' and 'DevType_y'. With Maps_Cartesian, we represent
% the velocity field as function of iy and ix indices. Look at the velocity
% field and give coherence with real field, then with the derivative.
%% Along the y axis, 
% the indices in the matrix are increasing downward on the vertical, if the 
% y velocity component is positive upward/downward then the derivative d/dy
%  must be descend/ascend along iy indices (it reverses axis and indices).
%% Along the x axis, 
% the indices in the matrix are increasing rightward on the horizontal, if  
% the x velocity component is positive rightward/leftward then the 
% derivative d/dy must be ascend/descend along iy indices (it reverses axis
% and idices).
%% => d/dx and d/dy must follow the direction of x and y axis given by the 
%%    positive Vx and Vy velocity component. 
%     1/ define the x and y axis: make your flow corresponds to its real
%     image, then positive velocity gives the direction of the axis.
%     2/ define the ascendent/descendent derivative along indices in order
%     to follow the x and y axis positive direction.
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
[returnOK] = Maps_Cartesian(Ux,Uy,XX,YY,idxmin,idxmax,idymin,idymax,DevType_x,DevType_y)
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