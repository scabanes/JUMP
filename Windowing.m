classdef Windowing
    methods(Static)
%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   HANNING WINDOWING
%  ###########################################################################################################################
%  ###########################################################################################################################
function H = WindowingHanning(Y,dim)
%  INPUT: N la taille du signal, Y le signal de N echantillon en temps ou espace.
% 	 L longueur du signal dans  l unite souhaitee.
%  OUTPUT: Y_Han de taille N 
% 
%  % % Y*H(n) % on va perdre ici en amplitude du signal, si AMP=1 dans H
%  ==> l'apmplitude est comprise entre 0 (au bord) et 1.
%  avec Hanning on va diminuer l amplitude du signal par 2, integral de la
%  courbe une fois passee dans la fenetre de hanning en revanche le pic en
%  frq est plus marque, le signal etant periodisise.
%  ==> En revanche les basses frquencs sont plus energetiques (hanning ajoute
%  artificiellement un signal basse frequence) mais les hautes frqce sont
%  mieux resolues.
%  ==> si le signal est deja periodique on gagne rien avec hanning on perd
%  au contraire les resolution basse frq en revanche si non periodique on
%  perd bcp dans une fft classique

%  ==> hanning function
    H=ones(size(Y));
if(dim==1)
    N=size(Y,1);
    Han=0.5*(1.-cos(2.*pi*[0:1:N-1]'/(N-1.)));
    H=H.*Han;
else
    N=size(Y,2);
    Han=0.5*(1.-cos(2.*pi*[0:1:N-1]'/(N-1.)));
    H=H.*Han';
end

end

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   TUKEY WINDOWING
%  ###########################################################################################################################
%  ###########################################################################################################################

function T = WindowingTukey(Y,dim)
% %  INPUT: N la taille du signal, Y le signal de N echantillon en temps ou espace.
% 	 L longueur du signal dans  l unite souhaitee.
%  OUTPUT: Y_Han de taille N 
% 
%  % % Y*H(n) % on va perdre ici en amplitude du signal, si AMP=1 dans H
%  ==> l'apmplitude est comprise entre 0 (au bord) et 1.
%  avec Hanning on va diminuer l amplitude du signal par 2, integral de la
%  courbe une fois passee dans la fenetre de hanning en revanche le pic en
%  frq est plus marque, le signal etant periodisise.
%  ==> En revanche les basses frquencs sont plus energetiques (hanning ajoute
%  artificiellement un signal basse frequence) mais les hautes frqce sont
%  mieux resolues.
%  ==> si le signal est deja periodique on gagne rien avec hanning on perd
%  au contraire les resolution basse frq en revanche si non periodique on
%  perd bcp dans une fft classique
T=ones(size(Y));
if(dim==1)
    N=size(Y,1);
%  ==> Tukey function
	Tukey = zeros(N,1);
	alpha = 0.05;
	n = [0:1:N-1]';
	n1 = alpha*(N-1)/2.;
	n2 = (N-1)*(1-alpha/2.);
	for iN=1:N
		if(n(iN) < n1)
			Tukey(iN) = 0.5*(1.+cos(pi*(2*n(iN)/(alpha*(N-1.)) - 1.)));
        elseif(n(iN) <= n2)
			Tukey(iN) = 1.;
        else
			Tukey(iN) = 0.5*(1.+cos(pi*(2*n(iN)/(alpha*(N-1.)) - 2./alpha + 1.)));
        end
    end
T=T.*Tukey;
else
    N=size(Y,2);
%  ==> Tukey function
	Tukey = zeros(N,1);
	alpha1 = 0.05;
    alpha2 = 0.05;
	n = [0:1:N-1]';
	n1 = alpha1*(N-1)/2.;
	n2 = (N-1)*(1-alpha2/2.);
	for iN=1:N
		if(n(iN) < n1)
			Tukey(iN) = 0.5*(1.+cos(pi*(2*n(iN)/(alpha1*(N-1.)) - 1.)));
        elseif(n(iN) <= n2)
			Tukey(iN) = 1.;
        else
			Tukey(iN) = 0.5*(1.+cos(pi*(2*n(iN)/(alpha2*(N-1.)) - 2./alpha2 + 1.)));
        end
    end
T=T.*Tukey';
end

end
    end
end
