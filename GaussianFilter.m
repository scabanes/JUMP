function [Field_fil_time_l,t_max] = GaussianFilter(Field,Vl,N_l,nodi,Tmax,Itime,passo_fr)
%____________________________________________________________________________________
% sig = l/sqrt(12). This relation is given considerating a gaussian function
% of the type g = 1/Norm exp[ -0.5 ((r-0)/sig)² ] with 
% Norm = Sum_r1 Sum_r2 G(r) where r1 and r2 are the two grid directions. One 
% may find that common low-pass filter [1-2] are G ~ exp[ -6 ((r-0)/l)² ] leading
% to the relation sig = l/sqrt(12). This definition give the scale l at
% which the gaussian is almost zero, which is much larger than the standard
% deviation sig. The gaussian expression for g is the one used by matlab
% functions. Then the filtred field is obtained by the convolution of
% conv(G,field). We note that l is a scale in grid points that goes up to
% nodi, the grid size. It has then to be converted in cm to obtain pysical
% scales. Good papers to read are [3-4].
%____________________________________________________________________________________
% REFERENCES:
% [1] Physical mechanism of the inverse energy cascade of two-dimensional turbulence: 
%     a numerical investigation Z. XIAO et al (2009), doi:10.1017/S0022112008004266
% [2] Turbulence : the filtering approach, Germano (1991)
% [3] Altimetric measurements of rotating turbulence: cyclone-anticyclone asymmetry, 
%     inertial and Kelvin waves and spectral characteristics, Y. D. Afanasyev
%     (2013), [http://dx.doi.org/10.1063/1.4826477]
% [4] Evidence for the double cascade scenario in two-dimensional turbulence
%     G. Boffetta and S. Musacchio (2010), DOI: 10.1103/PhysRevE.82.016307
%____________________________________________________________________________________
t_max= (Itime+Tmax-1)/passo_fr;
Field_fil_time_l = zeros (nodi*nodi,t_max,N_l);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               GAUSSIAN:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for il = 1:N_l

    l=Vl(il);
% Here we calculate sigma
    sig = l/(12).^0.5;
   
%-->>>  This is designed to emulate a Gaussian function used by filter2  
% % %   if sig <= 1
% % %     
% % %     h = fspecial ('gaussian',5,sig);
% % %    
% % %   else
% % %     
% % %     dim_max = round(2*(3*sig));
% % %     
% % %     if mod(dim_max,2)  ==  0
% % %   
% % %        %disp('dimensione pari');
% % %    
% % %          W = dim_max + 1;
% % %         
% % %     h = fspecial ('gaussian',W,sig);
% % %     
% % %     else
% % %         
% % %      W = dim_max  ;
% % %      h = fspecial ('gaussian',W,sig);
% % %    
% % %     end
% % %     
% % %      
% % %   end
  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              FILTERING:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
appo=0;
for t = Itime : passo_fr : Itime+Tmax-1
    
    appo=appo+1;
% We reshape field for each time
    Field_i =reshape(Field(:,t),nodi,nodi);
% The field is filtrated with the gaussian h
%   Field_fil = filter2(h,Field_i);
    Field_fil = imgaussfilt(Field_i,sig);
% The field is saved for each time
    Field_fil_time_l (:,appo,il) = Field_fil (:);   

end

end

disp('filtraggio ok')
end

